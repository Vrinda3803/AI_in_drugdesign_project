const searchInput = document.getElementById('search');
const resultsBox = document.getElementById('results');
const recordView = document.getElementById('record-view');
const loadState = document.getElementById('load-state');
let records = [];
let selected = null;
let generatedCandidateIndex = 0;

const esc = (value) => String(value ?? '').replace(/[&<>"']/g, (char) => ({
  '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'
}[char]));

function renderRecord(record, generated = false) {
  selected = record;
  recordView.hidden = false;
  const outputBanner = document.getElementById('generated-output-banner');
  outputBanner.hidden = !generated;
  if (generated) {
    document.getElementById('generated-method').textContent = 'Source-informed analogue based on published B2RP substitutions.';
    document.getElementById('generation-status').textContent = `Showing ${record.name}; sequence and evidence are loaded below.`;
  }
  document.getElementById('source-label').textContent = record.source.toUpperCase();
  document.getElementById('record-id').textContent = record.record_id ? `· ${record.record_id}` : '';
  document.getElementById('peptide-name').textContent = record.name;
  document.getElementById('peptide-aliases').textContent = record.aliases?.length
    ? `Also listed as ${record.aliases.join(' · ')}` : '';

  const candidateBadge = document.getElementById('candidate-badge');
  candidateBadge.hidden = record.source !== 'Source-informed hypothesis';
  document.getElementById('sequence').textContent = record.sequence;
  document.getElementById('length-label').textContent = `${record.sequence.length} residues`;
  const mods = record.modification_note || {};
  const notes = [];
  if (mods.n_terminus) notes.push(`N-terminus: ${mods.n_terminus}`);
  if (mods.c_terminus) notes.push(`C-terminus: ${mods.c_terminus}`);
  if (mods.unusual_amino_acids?.length) notes.push('Nonstandard residues present');
  document.getElementById('modification-note').textContent = notes.join(' · ');

  const score = record.docking_score;
  const clusterDensity = record.docking_cluster_density;
  document.getElementById('docking-status').textContent = record.docking_status || 'Not run';
  document.getElementById('docking-score-label').textContent = clusterDensity ? 'TOP-CLUSTER DENSITY' : 'DOCKING SCORE';
  document.getElementById('docking-score').textContent = clusterDensity
    ? String(clusterDensity)
    : score === null || score === undefined || score === '' ? 'Unavailable' : String(score);
  document.getElementById('docking-note').textContent = clusterDensity || score === null || score === undefined || score === ''
    ? (record.docking_note || 'Peptide docking has not been run. No score is available for this record.')
    : `Reported value for ${record.docking_target || 'the listed target'}.`;
  const dockingJobLink = document.getElementById('docking-job-link');
  dockingJobLink.hidden = !record.docking_job_url;
  if (record.docking_job_url) dockingJobLink.href = record.docking_job_url;

  renderToxicity(record.toxicity_observations || [], record.toxicity_predictions || []);
  renderActivity(record.activity_observations || []);
  resultsBox.hidden = true;
}

function renderActivity(observations) {
  const empty = document.getElementById('activity-empty');
  const summary = document.getElementById('activity-summary');
  const tableWrap = document.getElementById('activity-table-wrap');
  const caveat = document.getElementById('activity-caveat');
  document.getElementById('activity-count').textContent = `${observations.length} ${observations.length === 1 ? 'record' : 'records'}`;
  empty.hidden = observations.length > 0;
  summary.hidden = observations.length === 0;
  tableWrap.hidden = observations.length === 0;
  caveat.hidden = observations.length === 0;
  if (!observations.length) { tableWrap.replaceChildren(); return; }
  summary.innerHTML = `<strong>${observations.length} source-reported MIC records</strong> · Acinetobacter baumannii`;
  const ordered = [...observations].sort((a, b) => a.target.localeCompare(b.target) || a.pubmed_id.localeCompare(b.pubmed_id));
  const table = (rows) => `<div class="table-wrap"><table><thead><tr><th>ISOLATE / TARGET</th><th>REPORTED RESULT</th><th>CONCENTRATION</th><th>MEDIUM</th><th>INOCULUM</th><th>SOURCE</th></tr></thead><tbody>${rows.map((row) => {
    const paper = row.pubmed_id
      ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${esc(row.pubmed_id)}/" target="_blank" rel="noreferrer">PMID ${esc(row.pubmed_id)}</a>`
      : '<span>Unresolved</span>';
    return `<tr><td>${esc(row.target)}</td><td class="result-cell">${esc(row.reported_result)}</td><td>${esc(row.concentration)} <span class="unit-cell">${esc(row.unit)}</span></td><td>${esc(row.medium)}</td><td>${esc(row.cfu)}</td><td class="paper-cell">${paper}</td></tr>`;
  }).join('')}</tbody></table></div>`;
  const visible = ordered.slice(0, 6);
  const more = ordered.length > visible.length
    ? `<details class="all-measurements"><summary>View all ${ordered.length} MIC records</summary>${table(ordered.slice(6))}</details>` : '';
  tableWrap.innerHTML = table(visible) + more;
}

function renderToxicity(observations, predictions) {
  const predictionEmpty = document.getElementById('tox-prediction-empty');
  const predictionWrap = document.getElementById('tox-predictions-wrap');
  predictionEmpty.hidden = predictions.length > 0;
  predictionWrap.hidden = predictions.length === 0;
  if (predictions.length) {
    predictionWrap.innerHTML = `<table><thead><tr><th>TOOL</th><th>ENDPOINT</th><th>PREDICTION</th><th>REPORTED SCORE</th><th>PPV</th><th>MODEL DETAILS</th></tr></thead><tbody>${predictions.map((row) => {
      const score = row.score ? `${esc(row.score)} ${esc(row.score_unit)}`.trim() : 'Not reported';
      const details = [row.model, row.threshold ? `threshold ${row.threshold}` : '', row.run_date].filter(Boolean).map(esc).join(' · ');
      const scoreDisplay = row.score_name ? `${esc(row.score_name)}: ${score}` : score;
      return `<tr><td>${esc(row.predictor)}</td><td>${esc(row.endpoint)}</td><td class="result-cell">${esc(row.prediction)}</td><td>${scoreDisplay}</td><td>${row.ppv ? esc(row.ppv) : '—'}</td><td>${details || 'Not specified'}</td></tr>`;
    }).join('')}</tbody></table>`;
  } else {
    predictionWrap.replaceChildren();
  }

  const empty = document.getElementById('tox-empty');
  const summary = document.getElementById('tox-summary');
  const tableWrap = document.getElementById('tox-table-wrap');
  const caveat = document.getElementById('tox-caveat');
  document.getElementById('tox-count').textContent = `${observations.length} ${observations.length === 1 ? 'record' : 'records'}`;
  empty.hidden = observations.length > 0;
  summary.hidden = observations.length === 0;
  tableWrap.hidden = observations.length === 0;
  caveat.hidden = observations.length === 0;
  if (!observations.length) return;

  const rbc = observations.filter((row) => /erythrocyte|red blood cell/i.test(row.target_cell));
  const cyto = observations.filter((row) => /cytotox|cell death|killing/i.test(row.endpoint_group));
  summary.innerHTML = `<strong>${observations.length} source-reported assays</strong> · ${rbc.length} erythrocyte/hemolysis records · ${cyto.length} cytotoxicity/cell-effect records`;

  const ordered = [...observations].sort((a, b) => {
    const rank = (x) => /human erythrocyte/i.test(x.target_cell) && /hemolysis/i.test(x.endpoint_group) ? 0
      : /hemolysis/i.test(x.endpoint_group) ? 1 : /cytotox|cell death/i.test(x.endpoint_group) ? 2 : 3;
    return rank(a) - rank(b) || a.target_cell.localeCompare(b.target_cell);
  });
  const visible = ordered.slice(0, 5);
  const table = (rows) => `<div class="table-wrap"><table><thead><tr><th>TARGET CELL</th><th>ASSAY</th><th>REPORTED RESULT</th><th>DOSE</th><th>SOURCE</th></tr></thead><tbody>${rows.map((row) => {
    const paper = row.pubmed_id
      ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${esc(row.pubmed_id)}/" target="_blank" rel="noreferrer">PMID ${esc(row.pubmed_id)}</a>`
      : '<span>Unresolved</span>';
    return `<tr><td>${esc(row.target_cell)}</td><td>${esc(row.endpoint_group)}</td><td class="result-cell">${esc(row.reported_result)}</td><td>${esc(row.concentration)} <span class="unit-cell">${esc(row.unit)}</span></td><td class="paper-cell">${paper}</td></tr>`;
  }).join('')}</tbody></table></div>`;
  const more = ordered.length > visible.length
    ? `<details class="all-measurements"><summary>View all ${ordered.length} measurements</summary>${table(ordered.slice(5))}</details>` : '';
  tableWrap.innerHTML = table(visible) + more;
}

function search(query) {
  const q = query.trim().toLocaleLowerCase();
  if (!q) { resultsBox.hidden = true; return; }
  const matches = records.map((record) => {
    const values = [record.name, ...(record.aliases || []), record.sequence, record.id, record.record_id, record.source].map((x) => String(x || '').toLocaleLowerCase());
    let score = 99;
    if (values.some((value) => value === q)) score = 0;
    else if (values.some((value) => value.startsWith(q))) score = 1;
    else if (values.some((value) => value.includes(q))) score = 2;
    return { record, score };
  }).filter((item) => item.score < 99)
    .sort((a, b) => a.score - b.score || a.record.name.localeCompare(b.record.name))
    .slice(0, 12);

  resultsBox.replaceChildren();
  if (!matches.length) {
    const empty = document.createElement('div');
    empty.className = 'no-results';
    empty.textContent = 'No matching peptide. Try a name, sequence fragment, or record ID.';
    resultsBox.append(empty);
  } else {
    for (const { record } of matches) {
      const button = document.createElement('button');
      button.type = 'button';
      button.className = 'result-option';
      button.setAttribute('role', 'option');
      button.innerHTML = `<span><span class="result-name">${esc(record.name)}</span><br><span class="result-seq">${esc(record.sequence)}</span></span><span class="result-source">${esc(record.source)} · ${esc(record.record_id)}</span>`;
      button.addEventListener('click', () => {
        renderRecord(record);
        searchInput.value = record.name;
      });
      resultsBox.append(button);
    }
  }
  resultsBox.hidden = false;
}

searchInput.addEventListener('input', () => search(searchInput.value));
document.getElementById('clear-search').addEventListener('click', () => {
  searchInput.value = '';
  resultsBox.hidden = true;
  searchInput.focus();
});
document.getElementById('generate-candidate').addEventListener('click', () => {
  const bacterium = document.getElementById('bacterium-select').value;
  const candidates = records.filter((record) => record.candidate_status === 'generated_unvalidated_hypothesis');
  if (bacterium !== 'Acinetobacter baumannii' || !candidates.length) {
    document.getElementById('generation-status').textContent = 'No candidate generation panel is available for this selection.';
    return;
  }
  const candidate = candidates[generatedCandidateIndex % candidates.length];
  generatedCandidateIndex += 1;
  renderRecord(candidate, true);
  searchInput.value = candidate.name;
  recordView.scrollIntoView({ behavior: 'smooth', block: 'start' });
});
document.getElementById('copy-sequence').addEventListener('click', async (event) => {
  if (!selected) return;
  const button = event.currentTarget;
  try {
    await navigator.clipboard.writeText(selected.sequence);
    button.innerHTML = '<span aria-hidden="true">✓</span> Copied';
    setTimeout(() => { button.innerHTML = '<span aria-hidden="true">▣</span> Copy'; }, 1300);
  } catch {
    button.textContent = 'Select sequence to copy';
  }
});
document.addEventListener('click', (event) => {
  if (!event.target.closest('.lookup-panel')) resultsBox.hidden = true;
});

const data = window.PEPTIDE_DATA;
if (data?.records?.length) {
  records = data.records;
  const snapshotDate = new Date(`${data.snapshot_date}T00:00:00Z`);
  document.getElementById('snapshot-date').textContent = Number.isNaN(snapshotDate.getTime())
    ? data.snapshot_date
    : new Intl.DateTimeFormat('en-GB', { day: '2-digit', month: 'short', year: 'numeric', timeZone: 'UTC' })
      .format(snapshotDate).toUpperCase();
  loadState.hidden = true;
  const preferred = records.find((record) => record.id === 'DBAASP:705')
    || records.find((record) => record.id === 'Generated:D4R_hypothesis')
    || records[0];
  renderRecord(preferred);
  searchInput.value = preferred.name;
} else {
  loadState.textContent = 'Peptide data is missing. Rebuild the local snapshot with ui/build_ui_data.py.';
}
