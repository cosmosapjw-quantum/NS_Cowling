# Bounded post-change verification

The local full suites ran on clean implementation commit 594552d7a597dd150be14b9a860a8f5eaa571d7f.
The subsequent tracked test change only sets abs=0.0 on the near-full-volume assertion.
Its targeted suite was then rerun. The helper source is byte-identical throughout.

The complete Decimal comparison is retained compressed, including all 1,722 scalar outputs.
This is a new reproducible check, not recovery of the previously unretained 630-case sweep.
The source HTML and AME2016 identities in the intake report were not independently reacquired here.
The original nine logs remain unchanged and are covered by the parent SHA256SUMS.txt.

Replay from repository root into a new build directory:

~~~sh
mkdir -p build/review_scurto_replay
cp data/receipts/core_crust_followup_20261116/verification/decimal_check.py build/review_scurto_replay/check.py
PYTHONPATH=python_frontend python3 build/review_scurto_replay/check.py
~~~

The script writes its result alongside the copied script; keep the retained receipt immutable.
Replays record their actual HEAD. Source hash and output rows determine comparison identity.
To verify retained files, run sha256sum -c SHA256SUMS.txt from this directory.
CI #614 additionally reran the existing proxy census; this did not calibrate the helper.
No new stellar calculation was run. Scientific status remains PARTIALLY_COMPLETE.
