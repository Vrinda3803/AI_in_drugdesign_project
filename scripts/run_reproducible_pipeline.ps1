<#
.SYNOPSIS
  Rebuild the local AMP research artifacts from cached source data.
.DESCRIPTION
  Does not submit sequences to external predictors or start CABSdock. Use
  -RefreshDbaasp to fetch any missing DBAASP detail records, and -BuildUi to
  refresh the lookup snapshot after the analysis stages.
#>
[CmdletBinding()]
param(
    [switch]$RefreshDbaasp,
    [switch]$BuildUi,
    [switch]$SkipPyMol
)

$ErrorActionPreference = 'Stop'
$ProjectRoot = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $ProjectRoot
$StartedUtc = [DateTime]::UtcNow.ToString('o')
$CompletedStages = [System.Collections.Generic.List[string]]::new()

$Python = Join-Path $ProjectRoot '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $Python)) {
    throw "Project Python environment not found at $Python. Create it and install requirements.txt first."
}

function Invoke-ProjectScript {
    param([Parameter(Mandatory)][string]$RelativePath)
    $ScriptPath = Join-Path $ProjectRoot $RelativePath
    if (-not (Test-Path -LiteralPath $ScriptPath)) {
        throw "Required pipeline script is missing: $RelativePath"
    }
    Write-Host "`n>>> $RelativePath" -ForegroundColor Cyan
    & $Python $ScriptPath
    if ($LASTEXITCODE -ne 0) {
        throw "Pipeline stopped: $RelativePath exited with code $LASTEXITCODE"
    }
    $CompletedStages.Add($RelativePath)
}

if ($RefreshDbaasp) {
    Invoke-ProjectScript 'scripts/data_collection/download_dbaasp_ab_details.py'
}

$Stages = @(
    'scripts/data_collection/prepare_amp_datasets.py',
    'scripts/preprocessing/parse_mic_observations.py',
    'scripts/features/engineer_sequence_features.py',
    'scripts/preprocessing/extract_resistant_subcohort.py',
    'scripts/preprocessing/extract_paper_review_panel.py',
    'scripts/preprocessing/extract_toxicity_records.py',
    'scripts/preprocessing/prepare_human_rbc_h50_candidates.py',
    'scripts/preprocessing/build_activity_toxicity_evidence_panel.py',
    'scripts/activity_prediction/build_b2rp_lead_benchmark.py',
    'scripts/sequence_generation/generate_b2rp_hypotheses.py',
    'scripts/sequence_generation/export_candidate_fasta.py',
    'scripts/activity_prediction/build_multiobjective_candidate_priority.py',
    'scripts/activity_prediction/mic_regression_baselines.py',
    'scripts/activity_prediction/mic_relative_classifier.py',
    'scripts/activity_prediction/h50_regression_baselines.py',
    'scripts/docking/acquire_ftsZ_structures.py',
    'scripts/docking/summarize_ftsZ_inventory.py'
)
foreach ($Stage in $Stages) { Invoke-ProjectScript $Stage }

if (-not $SkipPyMol) {
    Invoke-ProjectScript 'scripts/docking/prepare_ftsZ_with_pymol.py'
} else {
    Write-Host 'Skipping PyMOL structure preparation by request.' -ForegroundColor Yellow
}

Invoke-ProjectScript 'scripts/docking/prepare_cabsdock_run_matrix.py'
Invoke-ProjectScript 'scripts/docking/summarize_cabsdock_runs.py'

if ($BuildUi) {
    Invoke-ProjectScript 'ui/build_ui_data.py'
}

$SafetyPredictionsPath = Join-Path $ProjectRoot 'results\candidate_prioritization\external_safety_predictions.csv'
$SafetyPredictionRows = if (Test-Path -LiteralPath $SafetyPredictionsPath) { @(Import-Csv -LiteralPath $SafetyPredictionsPath).Count } else { 0 }
$WebDockingJobsPath = Join-Path $ProjectRoot 'results\docking\cabsdock\web_jobs.csv'
$WebDockingJobs = if (Test-Path -LiteralPath $WebDockingJobsPath) { @(Import-Csv -LiteralPath $WebDockingJobsPath) } else { @() }
$DockingStatusCounts = @{}
foreach ($Group in ($WebDockingJobs | Group-Object -Property status)) {
    $DockingStatusCounts[$Group.Name] = $Group.Count
}
$SafetyStatus = if ($SafetyPredictionRows -gt 0) { "imported_externally; $SafetyPredictionRows normalized rows; not submitted by this pipeline" } else { 'not_run_by_this_pipeline' }
$DockingStatus = if ($WebDockingJobs.Count -gt 0) { "$($WebDockingJobs.Count) external web jobs recorded; see web_jobs.csv; not submitted by this pipeline" } else { 'run_plan_prepared_only; no poses or scores generated' }

$RunManifest = [ordered]@{
    started_utc = $StartedUtc
    finished_utc = [DateTime]::UtcNow.ToString('o')
    refresh_dbaasp_requested = [bool]$RefreshDbaasp
    ui_build_requested = [bool]$BuildUi
    pymol_skipped = [bool]$SkipPyMol
    completed_stages = @($CompletedStages)
    external_safety_predictions = $SafetyStatus
    external_safety_prediction_rows = $SafetyPredictionRows
    cabsdock = $DockingStatus
    cabsdock_web_job_count = $WebDockingJobs.Count
    cabsdock_web_job_status_counts = $DockingStatusCounts
}
$ManifestPath = Join-Path $ProjectRoot 'results\pipeline_run_manifest.json'
$ManifestDirectory = Split-Path -Parent $ManifestPath
New-Item -ItemType Directory -Force -Path $ManifestDirectory | Out-Null
$RunManifest | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $ManifestPath -Encoding utf8

Write-Host "`nLocal rebuild finished. The pipeline does not submit sequences to external predictors or start docking; it preserves and reports existing external predictions and docking-job links without fabricating results." -ForegroundColor Green
