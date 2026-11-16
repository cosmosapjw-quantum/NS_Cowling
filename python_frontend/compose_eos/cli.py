from __future__ import annotations
import argparse,json
from pathlib import Path
from .api import read_compose,check_consistency,normalize_compose

def main(argv=None):
    parser=argparse.ArgumentParser(prog="compose-eos")
    sub=parser.add_subparsers(dest="command",required=True)
    r=sub.add_parser("read");r.add_argument("directory",type=Path);r.add_argument("--out",type=Path)
    v=sub.add_parser("verify");v.add_argument("directory",type=Path);v.add_argument("--json",type=Path)
    n=sub.add_parser("normalize");n.add_argument("directory",type=Path);n.add_argument("output",type=Path)
    sub.add_parser("generators",help="List generation representations and physical scope")
    g=sub.add_parser("generate",help="Generate one explicitly scoped model from a source checkout")
    from .generation_dispatch import MODEL_IDS,SPG_VARIANTS
    g.add_argument("model",choices=MODEL_IDS);g.add_argument("output",type=Path)
    g.add_argument("--repository-root",type=Path)
    g.add_argument("--points",type=int)
    g.add_argument("--n-min",type=float);g.add_argument("--n-max",type=float)
    g.add_argument("--variant",choices=SPG_VARIANTS,default="AME2408_candidate")
    g.add_argument("--densities",nargs="+",type=float)
    g.add_argument("--resolutions",nargs="+",type=int)
    g.add_argument("--workers",type=int,default=1)
    args=parser.parse_args(argv)
    if args.command=="generators":
        from .generation_dispatch import generator_catalog
        print(json.dumps(generator_catalog(),indent=2));return 0
    if args.command=="generate":
        from .generation_dispatch import run_generation
        try:
            return run_generation(args.model,args.output,root=args.repository_root,
                points=args.points,n_min=args.n_min,n_max=args.n_max,variant=args.variant,
                densities=args.densities,resolutions=args.resolutions,workers=args.workers)
        except (ValueError,OSError) as exc:
            parser.error(str(exc))
    if args.command=="read":
        data=read_compose(args.directory)
        result={"status":"ok","n":data.n,"composition_available":data.composition_available,
                "nadd":data.q_extra_count,"source_dir":data.source_dir}
        if args.out:data.write_internal_directory(args.out)
        print(json.dumps(result,indent=2));return 0
    if args.command=="verify":
        data=read_compose(args.directory);report=check_consistency(data)
        text=json.dumps(report,indent=2,default=str)
        if args.json:args.json.write_text(text+"\n",encoding="utf-8")
        print(text);return 0 if report["overall"]["pass_"] else 1
    if args.command=="normalize":
        normalize_compose(args.directory,args.output);print(args.output);return 0
    return 2

if __name__=="__main__":raise SystemExit(main())
