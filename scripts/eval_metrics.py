#!/usr/bin/env python3
"""Ranking metrics for a folddisco result file, including average precision.

`folddisco benchmark` reports true positives at k false positives and recall;
it does not compute average precision, which is the headline metric in
feature_evaluation.md. This script computes both from the same inputs, so every
number in that document can be re-derived from the repository.

It replicates `benchmark --afdb-to-uniprot` exactly: an identifier is the
basename with its structure extension removed, then the second dash-separated
field (`AF-P17538-F1-model_v4.pdb` -> `P17538`), and duplicates are dropped
keeping first occurrence, so the AlphaFold fragments of one protein collapse to
a single accession. Ranking order is the order of lines in the result file.

Reported:
  hits_raw       lines in the result file (AFDB fragments, before deduplication)
  result_len     distinct accessions ranked
  tp@kfp         true positives found walking the ranked list until the k-th
                 false positive, for k in 1..500
  tp/fp/prec@topN at fixed depths N in 100..2000
  recall         tp_all / answer_len over the whole ranked list
  ap             average precision: for every rank i holding a true positive,
                 accumulate precision@i, then divide by the size of the FULL
                 answer set. Answers that never appear in the result therefore
                 contribute zero, which makes AP comparable between two
                 configurations that return different numbers of hits.

Usage:
  python3 scripts/eval_metrics.py <result.tsv> <answer.tsv> <index.lookup>

Reproduce a documented row (serine triad, --nonrigid, AP 0.8947):

  IDX=index/h_sapiens_folddisco
  folddisco query -i $IDX -p query/4CHA.pdb -q B57,B102,C195 -t 8 --nonrigid \
    --skip-match --per-structure --format-output tid > result.tsv
  python3 scripts/eval_metrics.py result.tsv data/serine_answer.tsv $IDX.lookup
"""
import sys, json

def parse(x):
    x = x.split('/')[-1]
    for ext in ('.pdb','.cif','.fcz','.ent'):
        if x.endswith(ext): x = x[:-4]; break
        if x.endswith(ext+'.gz'): x = x[:-7]; break
    s = x.split('-')
    return s[1] if len(s) >= 2 else x

def read_col(path, col=0):
    out, seen = [], set()
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line: continue
            v = parse(line.split('\t')[col])
            if v not in seen:
                seen.add(v); out.append(v)
    return out

def main(res, ans, lookup):
    result = read_col(res)
    answer = set(read_col(ans))
    allids = set(read_col(lookup, 1))
    n_ans = len(answer)
    o = {'hits_raw': sum(1 for _ in open(res)), 'result_len': len(result),
         'answer_len': n_ans, 'total_ids': len(allids)}
    # TP@kFP
    for k in (1,2,5,10,20,50,100,200,500):
        tp = fp = 0
        for t in result:
            if t in answer: tp += 1
            else: fp += 1
            if fp >= k: break
        o[f'tp@{k}fp'] = tp
    # precision@N / FP@N
    for n in (100,200,500,1000,2000):
        head = result[:n]
        tp = sum(1 for t in head if t in answer)
        o[f'tp@top{n}'] = tp
        o[f'fp@top{n}'] = len(head)-tp
        o[f'prec@top{n}'] = round(tp/len(head), 4) if head else None
    tp_all = sum(1 for t in result if t in answer)
    o['tp_all'] = tp_all; o['fp_all'] = len(result)-tp_all
    o['precision'] = round(tp_all/len(result),4) if result else None
    o['recall'] = round(tp_all/n_ans,4)
    # average precision over the ranked list (AP wrt full answer set)
    tp = 0; ap = 0.0
    for i,t in enumerate(result,1):
        if t in answer:
            tp += 1; ap += tp/i
    o['ap'] = round(ap/n_ans, 4)
    print(json.dumps(o))

if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
