"""Generate exactly four editable manuscript pages and their PDF.

No numerical result is populated from synthetic tests or incomplete scoring.
Images are placeholders until the primary scientific audit has completed.
ReportLab is used instead of requiring a TeX installation on Wahab.
"""
import json,html
from pathlib import Path
import pandas as pd
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor,white
from reportlab.platypus import Paragraph,Table,TableStyle,Image
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from common import ROOT,sha,load_protocol
W,H=612,792;M=48;CW=W-2*M;NAVY=HexColor('#153C4A');TEAL=HexColor('#007C83');GRAY=HexColor('#EEF3F5')
STYLE=ParagraphStyle('body',fontName='Times-Roman',fontSize=10.2,leading=13.5,spaceAfter=6,textColor=HexColor('#26353C'))
SMALL=ParagraphStyle('small',parent=STYLE,fontSize=8,leading=10.5)
HEAD=ParagraphStyle('head',parent=STYLE,fontName='Helvetica-Bold',fontSize=11,leading=14,textColor=NAVY)

def verified_status():
    p=ROOT/'results/audit_summary.json'
    if not p.exists():return None
    s=json.loads(p.read_text())
    if not s.get('real_model_scores') or s.get('AlphaGenome_in_primary_fit') or s.get('protocol_sha256')!=sha(ROOT/'config/protocol.json'):raise ValueError('Preprint must use the complete primary open-model audit only')
    for model in ['enformer','borzoi']:
        for shard in range(8):
            x=ROOT/'scores'/model/f'shard{shard:02d}.json'
            if not x.exists() or not json.loads(x.read_text()).get('complete'):raise ValueError('No complete native-model scoring proof')
    hashes=json.loads((ROOT/'results/result_CSV_sha256.json').read_text())
    for name,digest in hashes.items():
        if sha(ROOT/'results'/name)!=digest:raise ValueError('Result CSV changed: '+name)
    return s

def build():
    protocol=load_protocol();status=verified_status();out=ROOT/'preprint';out.mkdir(exist_ok=True)
    c=canvas.Canvas(str(out/'Workshop_Preprint.pdf'),pagesize=(W,H));c.setTitle('Deployment-time reliability audit of regulatory variant-effect models')
    pages=[];y=0
    def page(n,topic):
        nonlocal y
        if n>1:c.showPage()
        c.setFillColor(NAVY);c.rect(0,H-32,W,32,fill=1,stroke=0);c.setFillColor(white);c.setFont('Helvetica-Bold',9);c.drawString(M,H-21,'RELIABILITY AUDIT v2 | WORKSHOP DRAFT')
        c.setFillColor(NAVY);c.setFont('Helvetica',8);c.drawString(M,24,'Working draft. New model results '+('complete' if status else 'pending')+'. No clinical or causal claim.');c.drawRightString(W-M,24,f'{n} / 4')
        y=H-52;pages.append('\n\n<!-- PAGE '+str(n)+' -->\n## '+topic+'\n')
    def para(text,style=STYLE,gap=5):
        nonlocal y
        p=Paragraph(text,style);_,height=p.wrap(CW,1000)
        if y-height<45:raise ValueError('Manuscript page overflow; shorten text before exporting')
        p.drawOn(c,M,y-height);y-=height+gap;pages.append(text.replace('<b>','**').replace('</b>','**').replace('<br/>','\n')+'\n\n')
    def title(text):para(text,HEAD,7)
    def table(rows,widths):
        nonlocal y
        data=[[Paragraph(str(x),SMALL) for x in row] for row in rows]
        t=Table(data,colWidths=widths);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),GRAY),('LINEBELOW',(0,0),(-1,0),.6,NAVY),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('LINEBELOW',(0,1),(-1,-1),.2,HexColor('#CED9DE'))]))
        _,h=t.wrap(CW,1000);t.drawOn(c,M,y-h);y-=h+9
        pages.append('\n'.join(' | '.join(map(str,r)) for r in rows)+'\n\n')
    def figure(filename,label,caption,height=145):
        nonlocal y
        p=ROOT/'results'/filename
        if status and p.exists():
            im=Image(str(p));ratio=min(CW/im.imageWidth,height/im.imageHeight);im.drawWidth=im.imageWidth*ratio;im.drawHeight=im.imageHeight*ratio;im.drawOn(c,M+(CW-im.drawWidth)/2,y-im.drawHeight);y-=im.drawHeight+3
        else:
            c.setFillColor(GRAY);c.roundRect(M,y-height,CW,height,6,fill=1,stroke=0);c.setFillColor(NAVY);c.setFont('Helvetica-Bold',12);c.drawCentredString(W/2,y-height/2+10,label+' | RESULT PENDING');c.setFont('Helvetica',9);c.drawCentredString(W/2,y-height/2-8,'Autofilled only after complete native scoring and the locked audit.');y-=height+3
        para('<b>'+label+'.</b> '+caption,SMALL,7);pages.append('![ '+label+' ](../results/'+filename+')\n\n')
    page(1,'Question and evidence')
    para('Auditing deployment-time reliability of regulatory variant-effect models across independent MPRA studies',ParagraphStyle('title',parent=HEAD,fontSize=18,leading=21),9)
    para('Shadman Pathan | Working manuscript. Collaborators, authorship and venue remain to be confirmed.',SMALL,8)
    title('Abstract')
    abstract_tail='This draft contains actual public-data preparation and pending model-result panels; it reports no new predictive improvement.'
    if status:
        abstract_tail='The locked primary verdict was '+status['decision']+'. '
        for ds,z in status['datasets'].items():
            value=z['primary_error_reduction'];ci=z['CI97_5']
            if value is not None and all(v is not None for v in ci):abstract_tail+=f'{ds} sign-error reduction was {100*value:.2f} percentage points (97.5% CI {100*ci[0]:.2f} to {100*ci[1]:.2f}). '
    para('Variant-effect foundation models can correlate with reporter-assay measurements while making unreliable individual predictions. We ask whether model-output-only reliability scores improve selective prediction beyond effect magnitude on independent studies. The policy is fitted and calibrated on Kircher et al. saturation mutagenesis data, then frozen for Ulirsch K562 and Abell GM12878 MPRA validation. Open Enformer and Borzoi outputs support the primary analysis. Replicate-separated sign targets and approximate noise-ceiling correction distinguish inference reliability from assay dynamic range. Selection budgets are identical within genomic blocks; a conditional comparison also matches source-defined magnitude bins. An internally locked rule requires improvement on both validation studies. AlphaGenome fitting remains permission-gated and cannot determine the primary verdict. '+abstract_tail)
    title('1. Motivation and the negative evidence')
    para('The author-reported AlphaGenome pilot found that magnitude filtering improved prediction fidelity, while a cross-modality rule added only 0.004. After controlling for measured assay spread, prediction spread no longer predicted accuracy: partial rank correlation 0.19 (p=0.32) across 29 conditions and 0.09 across 21 elements. These prior results were not independently rerun here. They motivate a test of transferable reliability beyond magnitude and measurement noise, rather than a claim that prediction spread alone establishes reliability.')
    title('2. Public data and independent validation')
    qc=pd.read_csv(ROOT/'results/data_QC.csv').set_index('dataset')
    table([['Study / accession','Primary variants','Replication / role'],['Kircher 2019 / GSE126550',f"{int(qc.loc['Kircher','primary_variants']):,}",'21 elements, 19 grouped families; source fit/calibration; most variant SEs unavailable'],['Ulirsch 2016 / GSE87711',f"{int(qc.loc['Ulirsch','primary_variants']):,}",'K562 CTRL; 6 RNA and 2 DNA libraries; independent validation study'],['Abell 2022 / GSE174534',f"{int(qc.loc['Abell','primary_variants']):,}",'GM12878; 3 paired transfections; independent validation study']],[160,77,279])
    para('Counts, designs and allele mappings were anonymously downloaded and hash-verified. Validation excludes Kircher genomic blocks; Abell also excludes Ulirsch blocks. Repeated Abell constructs are canonicalized by source identifier before examining effects. MPRA-study independence does not establish independence from foundation-model pretraining.',SMALL)
    page(2,'Locked methods')
    title('3. Frozen model scores and deployment features')
    para('We use pinned, released Enformer and Borzoi PyTorch implementations and publicly downloadable CC BY 4.0 weights, with attribution. Enformer uses 196,608 bp input for its corrected TensorFlow-gamma conversion; Borzoi uses 524,288 bp. Each genomic REF/ALT contrast receives actual network inference in both orientations, in float32. Borzoi stranded tracks are explicitly remapped under reverse complement. Reference passes are reused only within fixed 2 kb context tiles. A 512 bp score window is scored near each variant; released target scaling is undone where identifiable. Borzoi output cropping is checked against its pointwise full head on Wahab.')
    para('The point forecast is the metadata-selected, cell-matched mean DNASE signed log2 fold score from Enformer. It is fixed before validation; no best-of-modality oracle is credited. Borzoi point accuracy is also reported. Candidate-block prediction spread, prediction magnitude, within-cell track dispersion, signed cross-model disagreement and missing-aware ATAC/CAGE/H3K27ac summaries are available at deployment. Measured assay spread, experimental effects, counts and replicates never enter the policy.')
    title('4. Source-only learning and calibration')
    para('Kircher canonical SNVs with at least 10 barcodes are grouped into 19 families, including shared GLOBIN and MYC groups. A fixed hash split assigns 14 families to fitting and 5 to calibration. A fixed histogram gradient-boosting classifier estimates sign-error risk on source variants with absolute published effect at least 0.5. Source-only isotonic calibration produces conditional error-probability estimates. Ranking uses raw classifier risk, so a monotonic calibration step is never credited with an incremental ranking gain. Hyperparameters, source scales and magnitude-quintile boundaries remain frozen across validation studies.')
    title('5. Variance-robust audit targets')
    para('ALT/REF effects are recomputed as log2 RNA allelic ratio minus log2 DNA allelic ratio with a fixed 0.5 pseudocount. For Abell, one paired replicate screens absolute effect at least 0.5 log2; the other two determine sign truth, rotating through three views. For Ulirsch, disjoint RNA halves use distinct DNA libraries, then swap. Each variant counts once with error averaged over qualifying views. Neither screen nor held-out truth affects selection. This reduces shared-noise selection bias but does not turn noisy replicate signs into certain biological truth.')
    para('Approximate mean-effect reliability: R = 1 - mean(SEM squared) / variance(mean effect). Corrected accuracy: r(model, mean effect) / sqrt(R). Abell SEM squared = variance(three paired effects) / 3. Ulirsch SEM squared = variance(six RNA ratios) / 6 + variance(two DNA ratios) / 2. Nonpositive R is unidentified; corrected correlations outside [-1,1] are unstable and excluded from adjusted-correlation audits, without clipping.',SMALL)
    title('6. Preregistered primary pass/fail rule')
    para('The plan is internally timestamped and hash-locked before v2 model outputs, not externally registered. Retain exactly floor(0.25 x n) per study/genomic block. PASS requires, in BOTH studies: at least 20 retained audited variants and 10 contributing blocks; positive error-reduction lower bound versus magnitude under a two-study 97.5% paired block-bootstrap CI; positive 95% lower bound within matched source-frozen magnitude bins; at least 90% of magnitude-baseline audited opportunities; and positive full-retention point-correlation 95% lower bound. Use 10,000 block resamples. Adequate evidence failing a benefit gate is FAIL; incomplete or underpowered evidence is INCONCLUSIVE. Secondary results cannot rescue FAIL.',SMALL)
    page(3,'Results')
    title('7. Primary external results')
    if status:
        para('<b>'+status['decision']+'.</b> '+('The locked open-model criterion passed on both external studies.' if status['decision']=='PASS' else 'The locked criterion did not establish reliable incremental improvement.' if status['decision']=='FAIL' else 'The locked criterion remains unevaluable or underpowered.'))
        comparisons=pd.read_csv(ROOT/'results/primary_bootstrap_comparisons.csv')
        table([['Study','Error reduction, pp','97.5% paired CI, pp']]+[[r.dataset,f'{100*r.error_reduction:.2f}',f'[{100*r.ci97_5_low:.2f}, {100*r.ci97_5_high:.2f}]'] for r in comparisons.itertuples()],[150,160,206])
    else:
        para('<b>PENDING.</b> Enformer/Borzoi inference has not run in the CPU-only authoring sandbox. Only public-data preparation, metadata contracts, native imports and software tests have completed. No mock prediction or fabricated number is a scientific result. New point accuracy, calibration and incremental reliability remain unknown.')
    figure('figure1_external_gain.png','Figure 1','Paired sign-error reduction versus magnitude at identical 25% retention on each independent study, with primary two-study 97.5% intervals. The CSV also reports 95% intervals. Positive values favor the source-trained policy.',160)
    figure('figure2_risk_retention.png','Figure 2','Replicate-separated sign accuracy across prespecified retention rates. Selected fractions are identical within blocks; audited effectful counts can differ and are reported explicitly.',170)
    para('Source-condition point accuracy, external per-study Pearson/Spearman correlation, classical noise-ceiling diagnostics, absolute precision with bootstrap intervals and source-calibrated probability reliability are exported as CSVs. Only the locked DNASE forecast and primary rule determine the main interpretation.',SMALL)
    page(4,'Confounding, limitations and references')
    title('8. Conditional reliability and assay dynamic range')
    figure('figure3_adjusted_reliability.png','Figure 3','Partial rank association of prediction spread and disagreement with block accuracy, controlling measured assay spread and replicate-derived mean SE. This is a fitted audit diagnostic, not an inference feature or proof of causal independence.',135)
    title('9. Interpretation and limitations')
    para('A positive result would support transferable selective ranking in these reporter assays, not clinical pathogenicity, calibrated uncertainty everywhere, or reliability of all foundation models. A negative result supports magnitude-only selection under the tested setting and argues against claiming unique reliability information from disagreement. Within-quintile matching controls magnitude coarsely, so residual continuous-magnitude effects remain possible. Absolute precision and calibration must accompany any ranking gain.')
    para('The primary claim concerns the fixed DNASE readout, not model-wide accuracy or a direct reporter-expression output. The studies differ in cell background, oligo length, promoter architecture and effect estimation. Source coefficients are not identical to reprocessed validation log ratios. Episomal reporter measurements do not recreate endogenous chromatin or full haplotypes. Models share training assays; loci may have appeared in pretraining. Source-family calibration is small and may not transfer. Prediction spread depends on the candidate set, so it is not a reliability claim for an isolated variant. Noise correction assumes additive unbiased errors; Ulirsch library covariance and common technical error are not fully identifiable. A 1 Mb block is a conservative clustering choice, not a proven independent LD unit. DNA/RNA sample quality rules and thresholds can change the target population; all exclusions remain auditable.')
    para('Open-model fitting uses only permissively released weights and outputs with attribution. AlphaGenome raw acquisition is optional. ANY normalization fitted on its outputs, risk training, probability calibration or partial-regression audit remains behind the explicit output-use permission gate. API access alone does not open that gate. Optional AlphaGenome analyses never alter the primary manuscript verdict.',SMALL)
    title('References and reproducibility')
    para('1. Kircher et al. (2019), Nature Communications. doi:10.1038/s41467-019-11526-w. GSE126550.<br/>2. Ulirsch et al. (2016), Cell. doi:10.1016/j.cell.2016.04.048. GSE87711.<br/>3. Abell et al. (2022), Science. doi:10.1126/science.abj5117. GSE174534.<br/>4. Avsec et al. (2021), Nature Methods. Enformer. doi:10.1038/s41592-021-01252-x.<br/>5. Linder et al. (2025), Nature Genetics. Borzoi. doi:10.1038/s41588-024-02053-6.<br/>6. Open ports: github.com/lucidrains/enformer-pytorch; github.com/johahi/borzoi-pytorch. Pinned commits, weight revisions, licenses, raw download URLs and checksums accompany the code.<br/>7. Seeds: 20261012 (models/split), 812071 (10,000 bootstrap resamples). Internal protocol SHA-256: '+sha(ROOT/'config/protocol.json')[:24]+'... Full digest in config/protocol_lock.json.',SMALL,0)
    c.save();(out/'Workshop_Preprint.md').write_text('# Deployment-time reliability audit\n'+''.join(pages))
    print('Wrote four-page draft; new inference status:',status['decision'] if status else 'PENDING')
if __name__=='__main__':build()
