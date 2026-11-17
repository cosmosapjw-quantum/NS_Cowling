"""C6 tests for the new source transport and explicit campaign scope."""
from pathlib import Path
import json
import shutil
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/"python_frontend"),str(ROOT/"scripts")]
from compose_eos.new_model_matching import (parse_provider,project_provider,
    check_origin_policy,NEW_PUBLIC)
from maxwell_three_family_campaign import check_same_family, valid_components
from run_new_model_maxwell_20261117 import make_pairs, NEW_NATIVE, OLD_CRUSTS

SOURCES=ROOT/"data/new_model_maxwell_20261117/sources"


class NewModelMatchingTests(unittest.TestCase):
    def test_raw_projection_preserves_all_seven_columns_and_auxiliary_tokens(self):
        for model in ("XMLSLZ_TM1","XMLSLZ_NL3"):
            with self.subTest(model=model), tempfile.TemporaryDirectory() as td:
                parsed=parse_provider(SOURCES/model)
                data,selection=project_provider(SOURCES/model,Path(td)/"out","crust")
                for j in range(7):
                    np.testing.assert_array_equal(getattr(data,f"Q{j+1}"),
                        parsed["q"][:parsed["core_start"],j])
                sidecar=json.loads(Path(selection["auxiliary_sidecar"]).read_text())
                self.assertEqual(sidecar["rows"][9]["auxiliary_tokens"],parsed["rows"][9][11:])
                self.assertEqual(sidecar["rows"][9]["raw_thermo_line"],parsed["raw_lines"][9])
                self.assertGreater(selection["muQ_vs_muE_max_abs_MeV"],0)
                self.assertEqual(selection["nadd_original"],5 if model=="XMLSLZ_TM1" else 7)
                self.assertTrue(data.n<parsed["n"].size)

    def test_rejects_wrong_nadd_and_phase_geometry(self):
        for damage in ("nadd","phase"):
            with self.subTest(damage=damage), tempfile.TemporaryDirectory() as td:
                target=Path(td)/"XMLSLZ_NL3";target.mkdir()
                for name in ("eos.nb","eos.thermo","eos.compo"):
                    shutil.copyfile(SOURCES/"XMLSLZ_NL3"/name,target/name)
                path=target/("eos.thermo" if damage=="nadd" else "eos.compo")
                lines=path.read_text().splitlines()
                i=1 if damage=="nadd" else 0
                tokens=lines[i].split()
                tokens[10 if damage=="nadd" else 3]="5" if damage=="nadd" else "2"
                lines[i]=" ".join(tokens);path.write_text("\n".join(lines)+"\n")
                with self.assertRaises(ValueError):parse_provider(target)

    def test_explicit_hybrid_and_historical_policy_remain_distinct(self):
        crust={"origin":"new_public","family":"new_public","model":"XMLSLZ_NL3","role":"crust"}
        core={"origin":"native_generated","family":"native_generated","model":"GM_GM1","role":"core"}
        self.assertEqual(check_origin_policy(crust,core),"reference_assisted_hybrid")
        with self.assertRaises(ValueError):check_same_family(crust,core)
        with self.assertRaises(ValueError):check_origin_policy({**crust,"model":"UNREGISTERED"},core)
        with self.assertRaises(ValueError):check_origin_policy(core,crust)

    def test_geometry_and_invalid_intervals_cannot_be_bridged(self):
        valid=[True,True,False,True,False,True,True]
        intervals=[{"left_row_zero_based":i,"right_row_zero_based":i+1,"valid":v} for i,v in enumerate(valid)]
        self.assertEqual(valid_components(intervals,8),[(0,3),(5,8)])

    def test_exact_bounded_120_unique_pairs(self):
        records=[]
        for origin,models,role,level in (
            ("native_generated",OLD_CRUSTS,"crust","fixed_old_source"),
            ("native_generated",("GM_GM1",),"core","fixed_old_source"),
            ("public_original",("GM_GM1",),"core","fixed_old_source"),
            ("new_native",NEW_NATIVE,"core","fine"),
            ("new_public",NEW_PUBLIC,"crust","fixed_public_precision"),
            ("new_public",NEW_PUBLIC[:4],"core","fixed_public_precision")):
            records += [{"id":f"{origin}:{m}:{role}:{level}","origin":origin,"model":m,"role":role,"resolution":level} for m in models]
        pairs=make_pairs(records)
        from collections import Counter
        self.assertEqual(Counter(x[0] for x in pairs),{
            "old_native_crust__new_native_core":65,"new_public_crust__old_public_GM1":5,
            "new_public_crust__new_public_core":20,"reference_assisted_hybrid":30})

if __name__=="__main__":unittest.main()
