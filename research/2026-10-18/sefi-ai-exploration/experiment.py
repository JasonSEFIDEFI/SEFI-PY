"""Structured-memory feasibility probe, not an AI or energy benchmark.

Inputs are already parsed (entity, attribute, value) events. Counts are logical
record accesses, NOT tokens, FLOPs, bytes, or joules. Both stores retain last
writer provenance. The ordinary indexed baseline deliberately has the same
semantics to expose which benefits are generic indexing benefits.
"""
from pathlib import Path
import random, json

def run(n, queries=1000):
    rng=random.Random(20261017+n)
    events=[(i%100, (i//100)%4, rng.randrange(1000000)) for i in range(n)]
    plain={}
    identity={}
    for version,(entity,attribute,value) in enumerate(events):
        plain[entity,attribute]=(value,version)
        # Proposed computational interpretation; not the physical SEFI equations.
        row=identity.setdefault(entity,{'origin':entity,'authorship':{},'sovereignty':entity,'warp':0})
        row['authorship'][attribute]=(value,version)
        row['warp']+=1
    keys=list(plain)
    scan_touches=0
    for _ in range(queries):
        key=rng.choice(keys)
        for version in range(n-1,-1,-1):
            scan_touches+=1
            entity,attribute,value=events[version]
            if (entity,attribute)==key:
                answer=(value,version)
                break
        assert answer==plain[key]==identity[key[0]]['authorship'][key[1]]
    # A last-value state cannot answer arbitrary historical queries.
    superseded=sum((value,i)!=plain[(e,a)] for i,(e,a,value) in enumerate(events))
    return {'events':n,'queries':queries,'entity_count':len(identity),'live_attributes':len(plain),
            'latest_query_agreement':1.0,'reverse_scan_record_reads':scan_touches,
            'indexed_value_reads':queries,'identity_value_reads':queries,
            'all_methods_ingestion_event_reads':n,
            'superseded_events_not_reconstructible_from_current_state':superseded,
            'history_events_if_preserved':n}

if __name__=='__main__':
    rows=[run(n) for n in (100,1000,10000,100000)]
    out=Path(__file__).parent
    (out/'results.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
    print(json.dumps(rows,indent=2))
