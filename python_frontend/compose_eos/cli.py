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
    args=parser.parse_args(argv)
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
