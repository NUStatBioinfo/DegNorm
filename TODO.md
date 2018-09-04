## DegNorm TODO

### Frank
- Test modules:
    - finish test_reads.py
    - nmf tests
    - gene processing tests
    - utils tests
    - visualization tests
- Check for pandoc install, attempt install, if pandoc available, convert HTML report to .PDF
- Ensure that genes unqualified for baseline selection are being correctly processed (e.g. scaled) - it's currently
possible for a gene's coverage curves to simply be scaled down if it's not sent through baseline selection at all.