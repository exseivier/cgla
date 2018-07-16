---
name: Bug in findOverlaps
about: findOverlaps returns NA values

---

**Describe the bug**
in get_intergenic_coordinates function it preforma a findOverlaps procedure of a flanked genomicRanges and mRNA sequences genomicRanges object and it returns NA values. I then removed them and I think that when I track back the names of the intergenic regions, gene names could be shuffled.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Nothing
