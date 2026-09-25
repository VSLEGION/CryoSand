"""Minimal runner for environments without pytest. Runs property-free tests only.
On a normal machine use `pytest -q` instead."""
import importlib
import sys
import traceback

sys.path.insert(0, ".")
failed = passed = 0
for modname in ["tests.test_limits"]:
    mod = importlib.import_module(modname)
    for name in sorted(n for n in dir(mod) if n.startswith("test_")):
        try:
            getattr(mod, name)()
            passed += 1
            print(f"PASS {name}")
        except Exception:  # noqa: BLE001
            failed += 1
            print(f"FAIL {name}")
            traceback.print_exc()
print(f"\n{passed} passed, {failed} failed")
sys.exit(1 if failed else 0)
