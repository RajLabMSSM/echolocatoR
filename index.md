# Developer

![](https://github.com/RajLabMSSM/echolocatoR/raw/master/inst/hex/hex.png "Hex sticker for echolocatoR")  
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
[![](https://img.shields.io/badge/doi-10.1093/bioinformatics/btab658-blue.svg)](https://doi.org/10.1093/bioinformatics/btab658)
[![](https://img.shields.io/badge/devel%20version-3.0.0-black.svg)](https://github.com/RajLabMSSM/echolocatoR)
[![](https://img.shields.io/github/languages/code-size/RajLabMSSM/echolocatoR.svg)](https://github.com/RajLabMSSM/echolocatoR)
[![](https://img.shields.io/github/last-commit/RajLabMSSM/echolocatoR.svg)](https://github.com/RajLabMSSM/echolocatoR/commits/master)  
[![R build
status](https://github.com/RajLabMSSM/echolocatoR/workflows/rworkflows/badge.svg)](https://github.com/RajLabMSSM/echolocatoR/actions)
[![](https://codecov.io/gh/RajLabMSSM/echolocatoR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/RajLabMSSM/echolocatoR)  
[![](https://codecov.io/gh/RajLabMSSM/echolocatoR/branch/master/graphs/icicle.svg "Codecov icicle graph")](https://app.codecov.io/gh/RajLabMSSM/echolocatoR/tree/master)  

#### Authors: *Brian Schilder, Jack Humphrey, Towfique Raj*

## `echolocatoR`: Automated statistical and functional fine-mapping with extensive access to genome-wide datasets.

### The *echoverse*

`echolocatoR` is part of the
[***echoverse***](https://github.com/topics/echoverse), a suite of R
packages designed to facilitate different steps in genetic fine-mapping.

`echolocatoR` calls each of these other packages (i.e. “modules”)
internally to create a unified pipeline. However, you can also use each
module independently to create your own custom workflows.

#### ***echoverse*** dependency graph

![echoverse dependency
graph](https://raw.githubusercontent.com/RajLabMSSM/echolocatoR/master/man/figures/echoverse.png)

> Made with [`echodeps`](https://github.com/RajLabMSSM/echodeps), yet
> another ***echoverse*** module. See [here for the interactive
> version](https://rajlabmssm.github.io/Fine_Mapping/echolocatoR.dep_graph.html)
> with package descriptions and links to each GitHub repo.

### *echoverse* status

| Package | Level | Version | CI Status | Code Coverage | Runners |
|:---|:---|:---|:---|:---|:---|
| [echolocatoR](https://github.com/RajLabMSSM/echolocatoR) | 4 | [![](https://img.shields.io/github/r-package/v/RajLabMSSM/echolocatoR?label=&color=black)](https://github.com/RajLabMSSM/echolocatoR) | [![rworkflows](https://github.com/RajLabMSSM/echolocatoR/actions/workflows/rworkflows.yml/badge.svg?branch=master)](https://github.com/RajLabMSSM/echolocatoR/actions/workflows/rworkflows.yml) | [![codecov](https://codecov.io/gh/RajLabMSSM/echolocatoR/branch/master/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/echolocatoR) | ![Ubuntu](https://img.shields.io/badge/-222222?style=flat-square&logo=ubuntu&logoColor=white)![macOS](https://img.shields.io/badge/-333333?style=flat-square&logo=apple&logoColor=white)![Windows](https://img.shields.io/badge/-444444?style=flat-square&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDMuNDQ5TDkuNzUgMi4xdjkuNDUxSDBtMTAuOTQ5LTkuNjAyTDI0IDB2MTEuNEgxMC45NDlNMCAxMi42aDkuNzV2OS40NTFMMCAyMC42OTlNMTAuOTQ5IDEyLjZIMjRWMjRsLTEyLjktMS44MDEiLz48L3N2Zz4=) |
| [catalogueR](https://github.com/RajLabMSSM/catalogueR) | 3 | [![](https://img.shields.io/github/r-package/v/RajLabMSSM/catalogueR?label=&color=black)](https://github.com/RajLabMSSM/catalogueR) | [![rworkflows](https://github.com/RajLabMSSM/catalogueR/actions/workflows/rworkflows.yml/badge.svg?branch=master)](https://github.com/RajLabMSSM/catalogueR/actions/workflows/rworkflows.yml) | [![codecov](https://codecov.io/gh/RajLabMSSM/catalogueR/branch/master/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/catalogueR) | ![Ubuntu](https://img.shields.io/badge/-222222?style=flat-square&logo=ubuntu&logoColor=white)![macOS](https://img.shields.io/badge/-333333?style=flat-square&logo=apple&logoColor=white)![Windows](https://img.shields.io/badge/-444444?style=flat-square&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDMuNDQ5TDkuNzUgMi4xdjkuNDUxSDBtMTAuOTQ5LTkuNjAyTDI0IDB2MTEuNEgxMC45NDlNMCAxMi42aDkuNzV2OS40NTFMMCAyMC42OTlNMTAuOTQ5IDEyLjZIMjRWMjRsLTEyLjktMS44MDEiLz48L3N2Zz4=) |
| [echoAI](https://github.com/RajLabMSSM/echoAI) | 3 | [![](https://img.shields.io/github/r-package/v/RajLabMSSM/echoAI?label=&color=black)](https://github.com/RajLabMSSM/echoAI) | [![rworkflows](https://github.com/RajLabMSSM/echoAI/actions/workflows/rworkflows.yml/badge.svg?branch=main)](https://github.com/RajLabMSSM/echoAI/actions/workflows/rworkflows.yml) | [![codecov](https://codecov.io/gh/RajLabMSSM/echoAI/branch/main/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/echoAI) | ![Ubuntu](https://img.shields.io/badge/-222222?style=flat-square&logo=ubuntu&logoColor=white)![macOS](https://img.shields.io/badge/-333333?style=flat-square&logo=apple&logoColor=white)![Windows](https://img.shields.io/badge/-444444?style=flat-square&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDMuNDQ5TDkuNzUgMi4xdjkuNDUxSDBtMTAuOTQ5LTkuNjAyTDI0IDB2MTEuNEgxMC45NDlNMCAxMi42aDkuNzV2OS40NTFMMCAyMC42OTlNMTAuOTQ5IDEyLjZIMjRWMjRsLTEyLjktMS44MDEiLz48L3N2Zz4=) |
| [echoannot](https://github.com/RajLabMSSM/echoannot) | 3 | [![](https://img.shields.io/github/r-package/v/RajLabMSSM/echoannot?label=&color=black)](https://github.com/RajLabMSSM/echoannot) | [![rworkflows](https://github.com/RajLabMSSM/echoannot/actions/workflows/rworkflows.yml/badge.svg?branch=main)](https://github.com/RajLabMSSM/echoannot/actions/workflows/rworkflows.yml) | [![codecov](https://codecov.io/gh/RajLabMSSM/echoannot/branch/main/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/echoannot) | ![Ubuntu](https://img.shields.io/badge/-222222?style=flat-square&logo=ubuntu&logoColor=white)![macOS](https://img.shields.io/badge/-333333?style=flat-square&logo=apple&logoColor=white)![Windows](https://img.shields.io/badge/-444444?style=flat-square&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDMuNDQ5TDkuNzUgMi4xdjkuNDUxSDBtMTAuOTQ5LTkuNjAyTDI0IDB2MTEuNEgxMC45NDlNMCAxMi42aDkuNzV2OS40NTFMMCAyMC42OTlNMTAuOTQ5IDEyLjZIMjRWMjRsLTEyLjktMS44MDEiLz48L3N2Zz4=) |
| [echofinemap](https://github.com/RajLabMSSM/echofinemap) | 3 | [![](https://img.shields.io/github/r-package/v/RajLabMSSM/echofinemap?label=&color=black)](https://github.com/RajLabMSSM/echofinemap) | [![rworkflows](https://github.com/RajLabMSSM/echofinemap/actions/workflows/rworkflows.yml/badge.svg?branch=master)](https://github.com/RajLabMSSM/echofinemap/actions/workflows/rworkflows.yml) | [![codecov](https://codecov.io/gh/RajLabMSSM/echofinemap/branch/master/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/echofinemap) | ![Ubuntu](https://img.shields.io/badge/-222222?style=flat-square&logo=ubuntu&logoColor=white)![macOS](https://img.shields.io/badge/-333333?style=flat-square&logo=apple&logoColor=white)![Windows](https://img.shields.io/badge/-444444?style=flat-square&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDMuNDQ5TDkuNzUgMi4xdjkuNDUxSDBtMTAuOTQ5LTkuNjAyTDI0IDB2MTEuNEgxMC45NDlNMCAxMi42aDkuNzV2OS40NTFMMCAyMC42OTlNMTAuOTQ5IDEyLjZIMjRWMjRsLTEyLjktMS44MDEiLz48L3N2Zz4=) |
| [echoplot](https://github.com/RajLabMSSM/echoplot) | 3 | [![](https://img.shields.io/github/r-package/v/RajLabMSSM/echoplot?label=&color=black)](https://github.com/RajLabMSSM/echoplot) | [![rworkflows](https://github.com/RajLabMSSM/echoplot/actions/workflows/rworkflows.yml/badge.svg?branch=master)](https://github.com/RajLabMSSM/echoplot/actions/workflows/rworkflows.yml) | [![codecov](https://codecov.io/gh/RajLabMSSM/echoplot/branch/master/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/echoplot) | ![Ubuntu](https://img.shields.io/badge/-222222?style=flat-square&logo=ubuntu&logoColor=white)![macOS](https://img.shields.io/badge/-333333?style=flat-square&logo=apple&logoColor=white)![Windows](https://img.shields.io/badge/-444444?style=flat-square&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDMuNDQ5TDkuNzUgMi4xdjkuNDUxSDBtMTAuOTQ5LTkuNjAyTDI0IDB2MTEuNEgxMC45NDlNMCAxMi42aDkuNzV2OS40NTFMMCAyMC42OTlNMTAuOTQ5IDEyLjZIMjRWMjRsLTEyLjktMS44MDEiLz48L3N2Zz4=) |
| [downloadR](https://github.com/RajLabMSSM/downloadR) | 2 | [![](https://img.shields.io/github/r-package/v/RajLabMSSM/downloadR?label=&color=black)](https://github.com/RajLabMSSM/downloadR) | [![rworkflows](https://github.com/RajLabMSSM/downloadR/actions/workflows/rworkflows.yml/badge.svg?branch=main)](https://github.com/RajLabMSSM/downloadR/actions/workflows/rworkflows.yml) | [![codecov](https://codecov.io/gh/RajLabMSSM/downloadR/branch/main/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/downloadR) | ![Ubuntu](https://img.shields.io/badge/-222222?style=flat-square&logo=ubuntu&logoColor=white)![macOS](https://img.shields.io/badge/-333333?style=flat-square&logo=apple&logoColor=white)![Windows](https://img.shields.io/badge/-444444?style=flat-square&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDMuNDQ5TDkuNzUgMi4xdjkuNDUxSDBtMTAuOTQ5LTkuNjAyTDI0IDB2MTEuNEgxMC45NDlNMCAxMi42aDkuNzV2OS40NTFMMCAyMC42OTlNMTAuOTQ5IDEyLjZIMjRWMjRsLTEyLjktMS44MDEiLz48L3N2Zz4=) |
| [echoconda](https://github.com/RajLabMSSM/echoconda) | 2 | [![](https://img.shields.io/github/r-package/v/RajLabMSSM/echoconda?label=&color=black)](https://github.com/RajLabMSSM/echoconda) | [![rworkflows](https://github.com/RajLabMSSM/echoconda/actions/workflows/rworkflows.yml/badge.svg?branch=main)](https://github.com/RajLabMSSM/echoconda/actions/workflows/rworkflows.yml) | [![codecov](https://codecov.io/gh/RajLabMSSM/echoconda/branch/main/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/echoconda) | ![Ubuntu](https://img.shields.io/badge/-222222?style=flat-square&logo=ubuntu&logoColor=white)![macOS](https://img.shields.io/badge/-333333?style=flat-square&logo=apple&logoColor=white) |
| [echoLD](https://github.com/RajLabMSSM/echoLD) | 2 | [![](https://img.shields.io/github/r-package/v/RajLabMSSM/echoLD?label=&color=black)](https://github.com/RajLabMSSM/echoLD) | [![rworkflows](https://github.com/RajLabMSSM/echoLD/actions/workflows/rworkflows.yml/badge.svg?branch=main)](https://github.com/RajLabMSSM/echoLD/actions/workflows/rworkflows.yml) | [![codecov](https://codecov.io/gh/RajLabMSSM/echoLD/branch/main/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/echoLD) | ![Ubuntu](https://img.shields.io/badge/-222222?style=flat-square&logo=ubuntu&logoColor=white)![macOS](https://img.shields.io/badge/-333333?style=flat-square&logo=apple&logoColor=white)![Windows](https://img.shields.io/badge/-444444?style=flat-square&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDMuNDQ5TDkuNzUgMi4xdjkuNDUxSDBtMTAuOTQ5LTkuNjAyTDI0IDB2MTEuNEgxMC45NDlNMCAxMi42aDkuNzV2OS40NTFMMCAyMC42OTlNMTAuOTQ5IDEyLjZIMjRWMjRsLTEyLjktMS44MDEiLz48L3N2Zz4=) |
| [echotabix](https://github.com/RajLabMSSM/echotabix) | 2 | [![](https://img.shields.io/github/r-package/v/RajLabMSSM/echotabix?label=&color=black)](https://github.com/RajLabMSSM/echotabix) | [![rworkflows](https://github.com/RajLabMSSM/echotabix/actions/workflows/rworkflows.yml/badge.svg?branch=main)](https://github.com/RajLabMSSM/echotabix/actions/workflows/rworkflows.yml) | [![codecov](https://codecov.io/gh/RajLabMSSM/echotabix/branch/main/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/echotabix) | ![Ubuntu](https://img.shields.io/badge/-222222?style=flat-square&logo=ubuntu&logoColor=white)![macOS](https://img.shields.io/badge/-333333?style=flat-square&logo=apple&logoColor=white) |
| [devoptera](https://github.com/RajLabMSSM/devoptera) | 1 | [![](https://img.shields.io/github/r-package/v/RajLabMSSM/devoptera?label=&color=black)](https://github.com/RajLabMSSM/devoptera) | [![rworkflows](https://github.com/RajLabMSSM/devoptera/actions/workflows/rworkflows.yml/badge.svg?branch=master)](https://github.com/RajLabMSSM/devoptera/actions/workflows/rworkflows.yml) | [![codecov](https://codecov.io/gh/RajLabMSSM/devoptera/branch/master/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/devoptera) | ![Ubuntu](https://img.shields.io/badge/-222222?style=flat-square&logo=ubuntu&logoColor=white)![macOS](https://img.shields.io/badge/-333333?style=flat-square&logo=apple&logoColor=white)![Windows](https://img.shields.io/badge/-444444?style=flat-square&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDMuNDQ5TDkuNzUgMi4xdjkuNDUxSDBtMTAuOTQ5LTkuNjAyTDI0IDB2MTEuNEgxMC45NDlNMCAxMi42aDkuNzV2OS40NTFMMCAyMC42OTlNMTAuOTQ5IDEyLjZIMjRWMjRsLTEyLjktMS44MDEiLz48L3N2Zz4=) |
| [echodata](https://github.com/RajLabMSSM/echodata) | 1 | [![](https://img.shields.io/github/r-package/v/RajLabMSSM/echodata?label=&color=black)](https://github.com/RajLabMSSM/echodata) | [![rworkflows](https://github.com/RajLabMSSM/echodata/actions/workflows/rworkflows.yml/badge.svg?branch=main)](https://github.com/RajLabMSSM/echodata/actions/workflows/rworkflows.yml) | [![codecov](https://codecov.io/gh/RajLabMSSM/echodata/branch/main/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/echodata) | ![Ubuntu](https://img.shields.io/badge/-222222?style=flat-square&logo=ubuntu&logoColor=white)![macOS](https://img.shields.io/badge/-333333?style=flat-square&logo=apple&logoColor=white)![Windows](https://img.shields.io/badge/-444444?style=flat-square&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDMuNDQ5TDkuNzUgMi4xdjkuNDUxSDBtMTAuOTQ5LTkuNjAyTDI0IDB2MTEuNEgxMC45NDlNMCAxMi42aDkuNzV2OS40NTFMMCAyMC42OTlNMTAuOTQ5IDEyLjZIMjRWMjRsLTEyLjktMS44MDEiLz48L3N2Zz4=) |
| [echogithub](https://github.com/RajLabMSSM/echogithub) | 1 | [![](https://img.shields.io/github/r-package/v/RajLabMSSM/echogithub?label=&color=black)](https://github.com/RajLabMSSM/echogithub) | [![rworkflows](https://github.com/RajLabMSSM/echogithub/actions/workflows/rworkflows.yml/badge.svg?branch=master)](https://github.com/RajLabMSSM/echogithub/actions/workflows/rworkflows.yml) | [![codecov](https://codecov.io/gh/RajLabMSSM/echogithub/branch/master/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/echogithub) | ![Ubuntu](https://img.shields.io/badge/-222222?style=flat-square&logo=ubuntu&logoColor=white)![macOS](https://img.shields.io/badge/-333333?style=flat-square&logo=apple&logoColor=white)![Windows](https://img.shields.io/badge/-444444?style=flat-square&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDMuNDQ5TDkuNzUgMi4xdjkuNDUxSDBtMTAuOTQ5LTkuNjAyTDI0IDB2MTEuNEgxMC45NDlNMCAxMi42aDkuNzV2OS40NTFMMCAyMC42OTlNMTAuOTQ5IDEyLjZIMjRWMjRsLTEyLjktMS44MDEiLz48L3N2Zz4=) |
| [echodeps](https://github.com/RajLabMSSM/echodeps) | \- | [![](https://img.shields.io/github/r-package/v/RajLabMSSM/echodeps?label=&color=black)](https://github.com/RajLabMSSM/echodeps) | [![rworkflows](https://github.com/RajLabMSSM/echodeps/actions/workflows/rworkflows.yml/badge.svg?branch=master)](https://github.com/RajLabMSSM/echodeps/actions/workflows/rworkflows.yml) | [![codecov](https://codecov.io/gh/RajLabMSSM/echodeps/branch/master/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/echodeps) | ![Ubuntu](https://img.shields.io/badge/-222222?style=flat-square&logo=ubuntu&logoColor=white)![macOS](https://img.shields.io/badge/-333333?style=flat-square&logo=apple&logoColor=white)![Windows](https://img.shields.io/badge/-444444?style=flat-square&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDMuNDQ5TDkuNzUgMi4xdjkuNDUxSDBtMTAuOTQ5LTkuNjAyTDI0IDB2MTEuNEgxMC45NDlNMCAxMi42aDkuNzV2OS40NTFMMCAyMC42OTlNMTAuOTQ5IDEyLjZIMjRWMjRsLTEyLjktMS44MDEiLz48L3N2Zz4=) |
| [echoverseTemplate](https://github.com/RajLabMSSM/echoverseTemplate) | \- | [![](https://img.shields.io/github/r-package/v/RajLabMSSM/echoverseTemplate?label=&color=black)](https://github.com/RajLabMSSM/echoverseTemplate) | [![rworkflows](https://github.com/RajLabMSSM/echoverseTemplate/actions/workflows/rworkflows.yml/badge.svg?branch=master)](https://github.com/RajLabMSSM/echoverseTemplate/actions/workflows/rworkflows.yml) | [![codecov](https://codecov.io/gh/RajLabMSSM/echoverseTemplate/branch/master/graph/badge.svg)](https://codecov.io/gh/RajLabMSSM/echoverseTemplate) | ![Ubuntu](https://img.shields.io/badge/-222222?style=flat-square&logo=ubuntu&logoColor=white)![macOS](https://img.shields.io/badge/-333333?style=flat-square&logo=apple&logoColor=white)![Windows](https://img.shields.io/badge/-444444?style=flat-square&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCI+PHBhdGggZmlsbD0id2hpdGUiIGQ9Ik0wIDMuNDQ5TDkuNzUgMi4xdjkuNDUxSDBtMTAuOTQ5LTkuNjAyTDI0IDB2MTEuNEgxMC45NDlNMCAxMi42aDkuNzV2OS40NTFMMCAyMC42OTlNMTAuOTQ5IDEyLjZIMjRWMjRsLTEyLjktMS44MDEiLz48L3N2Zz4=) |

> **Level** indicates dependency depth: **4** = top-level pipeline,
> **3** = analysis modules, **2** = mid-level utilities, **1** = base
> packages, **-** = non-pipeline. CI status reflects all runners
> combined. Click the badge to see per-runner results.

## Installation

#### Prerequisites

Before installing `echolocatoR`, make sure the following system
libraries are available. These are needed to compile R packages that
`echolocatoR` depends on.

**macOS** (via [Homebrew](https://brew.sh)):

``` bash
brew install libxml2 openssl curl zlib gcc
```

**Ubuntu / Debian**:

``` bash
sudo apt-get install -y libxml2-dev libssl-dev libcurl4-openssl-dev \
  zlib1g-dev gfortran
```

**Windows**: Install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R
version.

#### Install

``` r
if(!require("BiocManager")) install.packages("BiocManager")

BiocManager::install("RajLabMSSM/echolocatoR")
library(echolocatoR)
```

#### Diagnose issues

After installing, run the built-in diagnostic to check for any problems:

``` r
echolocatoR::check_echoverse_setup()
```

#### Troubleshooting

- **Dependency ordering**: Because `echolocatoR` relies on many
  subpackages that depend on one another, install errors can sometimes
  occur when R tries to update one package before its *echoverse*
  dependencies. The solution is to simply rerun
  `BiocManager::install("RajLabMSSM/echolocatoR")` until all subpackages
  are fully updated.
- **susieR version**: `echofinemap` requires `susieR` \>= 0.12.0, but
  CRAN sometimes has an older version. Fix:
  `devtools::install_github("stephenslab/susieR")`
- **GitHub rate limiting**: Installing many GitHub-hosted packages can
  trigger API rate limits. Set a GitHub Personal Access Token:
  `Sys.setenv(GITHUB_TOKEN = "your_token_here")`. See [GitHub
  docs](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)
  for how to create one.
- **Python / conda**: Some fine-mapping methods (PolyFun) and LD
  reference panels (UKB) require Python and conda. Install
  [Miniforge](https://github.com/conda-forge/miniforge) if you plan to
  use these features. The conda environment `echoR_mini` will be created
  automatically on first use.

### \[Optional\] [Docker/Singularity](https://rajlabmssm.github.io/echolocatoR/articles/docker)

`echolocatoR` now has its own dedicated Docker/Singularity container!
This greatly reduces issues related to system dependency conflicts and
provides a containerized interface for Rstudio through your web browser.
See [here for installation
instructions](https://rajlabmssm.github.io/echolocatoR/articles/docker).

## Quick start

``` r
library(echolocatoR)

## Load bundled fine-mapping results for the BST1 locus
dat <- echodata::BST1

## View consensus fine-mapped SNPs (agreed upon by multiple methods)
subset(dat, Consensus_SNP == TRUE,
       select = c(SNP, CHR, POS, P, mean.PP, Support))

## Visualize the locus
echoplot::plot_locus(
    dat = dat,
    locus_dir = file.path(tempdir(), echodata::locus_dir),
    LD_matrix = echodata::BST1_LD_matrix
)
```

See
[`vignette("echolocatoR")`](https://rajlabmssm.github.io/echolocatoR/articles/echolocatoR.md)
for the full pipeline, or
[`vignette("explore_results")`](https://rajlabmssm.github.io/echolocatoR/articles/explore_results.md)
to learn how to interpret results.

## Documentation

### [Website](https://rajlabmssm.github.io/echolocatoR)

### [Get started](https://rajlabmssm.github.io/echolocatoR/articles/echolocatoR)

### [Bugs/requests](https://github.com/RajLabMSSM/echolocatoR/issues)

Please report any bugs/requests on [GitHub
Issues](https://github.com/RajLabMSSM/echolocatoR/issues).

[Contributions](https://github.com/RajLabMSSM/echolocatoR/pulls) are
welcome!

### All *echoverse* vignettes

``` r
echoverse <- c('echolocatoR','echodata','echotabix',
               'echoannot','echoconda','echoLD',
               'echoplot','catalogueR','downloadR',
               'echofinemap','echodeps', # under construction
               'echogithub')
toc <- echogithub::github_pages_vignettes(owner = "RajLabMSSM",
                                          repo = echoverse,
                                          as_toc = TRUE,
                                          verbose = FALSE)
```

- ## 🦇 [echolocatoR](https://rajlabmssm.github.io/echolocatoR/)

  - ### [QTLs](https://rajlabmssm.github.io/echolocatoR//articles/QTLs.html)

  - ### [docker](https://rajlabmssm.github.io/echolocatoR//articles/docker.html)

  - ### [echolocatoR](https://rajlabmssm.github.io/echolocatoR//articles/echolocatoR.html)

  - ### [echoverse modules](https://rajlabmssm.github.io/echolocatoR//articles/echoverse_modules.html)

  - ### [explore results](https://rajlabmssm.github.io/echolocatoR//articles/explore_results.html)

  - ### [finemapping portal](https://rajlabmssm.github.io/echolocatoR//articles/finemapping_portal.html)

  - ### [plot locus](https://rajlabmssm.github.io/echolocatoR//articles/plot_locus.html)

  - ### [summarise](https://rajlabmssm.github.io/echolocatoR//articles/summarise.html)

- ## 🦇 [echodata](https://rajlabmssm.github.io/echodata/)

  - ### [echodata](https://rajlabmssm.github.io/echodata//articles/echodata.html)

  - ### [echolocatoR Finemapping Portal](https://rajlabmssm.github.io/echodata//articles/echolocatoR_Finemapping_Portal.html)

- ## 🦇 [echotabix](https://rajlabmssm.github.io/echotabix/)

  - ### [echotabix](https://rajlabmssm.github.io/echotabix//articles/echotabix.html)

- ## 🦇 [echoannot](https://rajlabmssm.github.io/echoannot/)

  - ### [cell type specific epigenomics](https://rajlabmssm.github.io/echoannot//articles/cell_type_specific_epigenomics.html)

  - ### [echoannot](https://rajlabmssm.github.io/echoannot//articles/echoannot.html)

- ## 🦇 [echoconda](https://rajlabmssm.github.io/echoconda/)

  - ### [echoconda](https://rajlabmssm.github.io/echoconda//articles/echoconda.html)

- ## 🦇 [echoLD](https://rajlabmssm.github.io/echoLD/)

  - ### [echoLD](https://rajlabmssm.github.io/echoLD//articles/echoLD.html)

- ## 🦇 [echoplot](https://rajlabmssm.github.io/echoplot/)

  - ### [echoplot](https://rajlabmssm.github.io/echoplot//articles/echoplot.html)

- ## 🦇 [catalogueR](https://rajlabmssm.github.io/catalogueR/)

  - ### [catalogueR](https://rajlabmssm.github.io/catalogueR//articles/catalogueR.html)

  - ### [colocalize](https://rajlabmssm.github.io/catalogueR//articles/colocalize.html)

- ## 🦇 [downloadR](https://rajlabmssm.github.io/downloadR/)

  - ### [downloadR](https://rajlabmssm.github.io/downloadR//articles/downloadR.html)

- ## 🦇 [echofinemap](https://rajlabmssm.github.io/echofinemap/)

  - ### [echofinemap](https://rajlabmssm.github.io/echofinemap//articles/echofinemap.html)

  - ### [echoverseTemplate](https://rajlabmssm.github.io/echofinemap//articles/echoverseTemplate.html)

- ## 🦇 [echodeps](https://rajlabmssm.github.io/echodeps/)

  - ### [echodeps](https://rajlabmssm.github.io/echodeps//articles/echodeps.html)

  - ### [echoverseTemplate](https://rajlabmssm.github.io/echodeps//articles/echoverseTemplate.html)

  - ### [my packages](https://rajlabmssm.github.io/echodeps//articles/my_packages.html)

- ## 🦇 [echogithub](https://rajlabmssm.github.io/echolocatoR/NA)

  - ### NA

## Introduction

Fine-mapping methods are a powerful means of identifying causal variants
underlying a given phenotype, but are underutilized due to the technical
challenges of implementation. `echolocatoR` is an R package that
automates end-to-end genomics fine-mapping, annotation, and plotting in
order to identify the most probable causal variants associated with a
given phenotype.

It requires minimal input from users (a GWAS or QTL summary statistics
file), and includes a suite of statistical and functional fine-mapping
tools. It also includes extensive access to datasets (linkage
disequilibrium panels, epigenomic and genome-wide annotations, QTL).

### Citation

If you use `echolocatoR`, or any of the **echoverse** modules, please
cite:

> Brian M Schilder, Jack Humphrey, Towfique Raj (2021) echolocatoR: an
> automated end-to-end statistical and functional genomic fine-mapping
> pipeline, *Bioinformatics*; btab658,
> <https://doi.org/10.1093/bioinformatics/btab658>

The elimination of data gathering and preprocessing steps enables rapid
fine-mapping of many loci in any phenotype, complete with locus-specific
publication-ready figure generation. All results are merged into a
single per-SNP summary file for additional downstream analysis and
results sharing. Therefore `echolocatoR` drastically reduces the
barriers to identifying causal variants by making the entire
fine-mapping pipeline rapid, robust and scalable.

![](https://raw.githubusercontent.com/RajLabMSSM/echolocatoR/master/man/figures/echolocatoR_Fig1.png)

## Literature

### For applications of `echolocatoR` in the literature, please see:

> 1.  E Navarro, E Udine, K de Paiva Lopes, M Parks, G Riboldi, BM
>     Schilder…T Raj (2020) Dysregulation of mitochondrial and
>     proteo-lysosomal genes in Parkinson’s disease myeloid cells.
>     Nature Genetics. <https://doi.org/10.1101/2020.07.20.212407>
> 2.  BM Schilder, T Raj (2021) Fine-Mapping of Parkinson’s Disease
>     Susceptibility Loci Identifies Putative Causal Variants. Human
>     Molecular Genetics, ddab294,
>     <https://doi.org/10.1093/hmg/ddab294>  
> 3.  K de Paiva Lopes, G JL Snijders, J Humphrey, A Allan, M Sneeboer,
>     E Navarro, BM Schilder…T Raj (2022) Genetic analysis of the human
>     microglial transcriptome across brain regions, aging and disease
>     pathologies. Nature Genetics,
>     <https://doi.org/10.1038/s41588-021-00976-y>

## `echolocatoR` v1.0 vs. v2.0

There have been a series of major updates between `echolocatoR` v1.0 and
v2.0. Here are some of the most notable ones (see **Details**):

- ***echoverse* subpackages**: `echolocatoR` has been broken into
  separate subpackages, making it much easier to edit/debug each step of
  the full `finemap_loci` pipeline, and improving robustness throughout.
  It also provides greater flexibility for users to construct their own
  custom pipelines from these modules.
- **`GITHUB_TOKEN`**: GitHub now requires users to create Personal
  Authentication Tokens (PAT) to avoid download limits. This is
  essential for installing `echolocatoR` as many resources from GitHub
  need to be downloaded. See
  [here](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)
  for further instructions. =
  [`echodata::construct_colmap()`](https://rdrr.io/pkg/echodata/man/construct_colmap.html):
  Previously, users were required to input key column name mappings as
  separate arguments to
  [`echolocatoR::finemap_loci`](https://rajlabmssm.github.io/echolocatoR/reference/finemap_loci.md).
  This functionality has been deprecated and replaced with a single
  argument, `colmap=`. This allows users to save the
  [`construct_colmap()`](https://rdrr.io/pkg/echodata/man/construct_colmap.html)
  output as a single variable and reuse it later without having to write
  out each mapping argument again (and helps reduce an already crowded
  list of arguments).
- **`MungeSumstats`**: `finemap_loci` now accepts the output of
  [`MungeSumstats::format_sumstats`/`import_sumstats`](https://github.com/neurogenomics/MungeSumstats)
  as-is (without requiring `colmap=`, so long as `munged=TRUE`).
  Standardizing your GWAS/QTL summary stats this way greatly reduces (or
  eliminates) the time taken to do manual formatting.
- **[`echolocatoR::finemap_loci`](https://rajlabmssm.github.io/echolocatoR/reference/finemap_loci.md)
  arguments**: Several arguments have been deprecated or had their names
  changed to be harmonized across all the subpackages and use a unified
  naming convention. See
  [`?echolocatoR::finemap_loci`](https://rajlabmssm.github.io/echolocatoR/reference/finemap_loci.md)
  for details.
- **`echoconda`**: The *echoverse* subpackage `echoconda` now handles
  all conda environment creation/use internally and automatically,
  without the need for users to create the conda environment themselves
  as a separate step. Also, the default conda env `echoR` has been
  replaced by `echoR_mini`, which reduces the number of dependencies to
  just the bare minimum (thus greatly speeding up build time and
  reducing potential version conflicts).
- **`FINEMAP`**: More outputs from the tool `FINEMAP` are now recorded
  in the `echolocatoR` results (see
  [`?echofinemap::FINEMAP`](https://rdrr.io/pkg/echofinemap/man/FINEMAP.html)
  or [this Issue](https://github.com/RajLabMSSM/echofinemap/issues/7)
  for details). Also, a common dependency conflict between
  `FINEMAP`\>=1.4 and MacOS has been resolved (see [this
  Issue](https://github.com/RajLabMSSM/echofinemap/issues/9) for
  details.
- **`echodata`**: All example data and data transformation functions
  have been moved to the *echoverse* subpackage
  [`echodata`](https://github.com/RajLabMSSM/echodata).
- **`LD_reference=`**: In addition to the *UKB*, *1KGphase1/3* LD
  reference panels,
  [`finemap_loci()`](https://rajlabmssm.github.io/echolocatoR/reference/finemap_loci.md)
  can now take custom LD panels by supplying
  `finemap_loci(LD_reference=)` with a list of paths to VCF files (.vcf
  / vcf.gz / vcf.bgz) or pre-computed LD matrices with RSIDs as the
  row/col names (.rda / .rds / .csv / .tsv. / .txt / .csv.gz / tsv.gz /
  txt.gz).
- **Expanded fine-mapping methods**: “ABF”, “COJO_conditional”,
  “COJO_joint” “COJO_stepwise”,“FINEMAP”,“PAINTOR” (including multi-GWAS
  and multi-ancestry fine-mapping),“POLYFUN_FINEMAP”
  ,“POLYFUN_SUSIE”,“SUSIE”
- **`FINEMAP` fixed**: There were a number of issues with `FINEMAP` due
  to differing output formats across different versions, system
  dependency conflicts, and the fact that it can produce multiple
  Credible Sets. All of these have been fixed and the latest version of
  `FINEMAP` can be run on all OS platforms.  
- **Debug mode**: Within
  [`finemap_loci()`](https://rajlabmssm.github.io/echolocatoR/reference/finemap_loci.md)
  I use a [`tryCatch()`](https://rdrr.io/r/base/conditions.html) when
  iterating across loci so that if one locus fails, the rest can
  continue. However this prevents using traceback feature in R, making
  debugging hard. Thus I now enabled debugging mode via a new argument:
  `use_tryCatch=FALSE`.

## Output descriptions

By default,
[`echolocatoR::finemap_loci()`](https://rajlabmssm.github.io/echolocatoR/reference/finemap_loci.md)
returns a nested list containing grouped by locus names (e.g. `$BST1`,
`$MEX3C`). The results of each locus contain the following elements:

- `finemap_dat`: Fine-mapping results from all selected methods merged
  with the original summary statistics (i.e. **Multi-finemap results**).
- `locus_plot`: A nested list containing one or more zoomed views of
  locus plots.  
- `LD_matrix`: The post-processed LD matrix used for fine-mapping.
- `LD_plot`: An LD plot (if used).
- `locus_dir`: Locus directory results are saved in.
- `arguments`: A record of the arguments supplied to `finemap_loci`.

In addition, the following object summarizes the results from the
locus-specific elements:  
- `merged_dat`: A merged `data.table` with all fine-mapping results from
all loci.

### Multi-finemap results files

The main output of `echolocatoR` are the multi-finemap files (for
example,
[`echodata::BST1`](https://rdrr.io/pkg/echodata/man/BST1.html)). They
are stored in the locus-specific *Multi-finemap* subfolders.

#### Column descriptions

- **Standardized GWAS/QTL summary statistics**: e.g.
  `SNP`,`CHR`,`POS`,`Effect`,`StdErr`. See `?finemap_loci()` for
  descriptions of each.  
- **leadSNP**: The designated proxy SNP per locus, which is the SNP with
  the smallest p-value by default.
- **\<tool\>.CS**: The 95% probability Credible Set (CS) to which a SNP
  belongs within a given fine-mapping tool’s results. If a SNP is not in
  any of the tool’s CS, it is assigned `NA` (or `0` for the purposes of
  plotting).  
- **\<tool\>.PP**: The posterior probability that a SNP is causal for a
  given GWAS/QTL trait.  
- **Support**: The total number of fine-mapping tools that include the
  SNP in its CS.
- **Consensus_SNP**: By default, defined as a SNP that is included in
  the CS of more than `N` fine-mapping tool(s), i.e. `Support>1`
  (default: `N=1`).  
- **mean.PP**: The mean SNP-wise PP across all fine-mapping tools used.
- **mean.CS**: If mean PP is greater than the 95% probability threshold
  (`mean.PP>0.95`) then `mean.CS` is 1, else 0. This tends to be a very
  stringent threshold as it requires a high degree of agreement between
  fine-mapping tools.

### Notes

- Separate multi-finemap files are generated for each LD reference panel
  used, which is included in the file name (e.g.
  *UKB_LD.Multi-finemap.tsv.gz*).
- Each fine-mapping tool defines its CS and PP slightly differently, so
  please refer to the associated original publications for the exact
  details of how these are calculated (links provided above).

## Fine-mapping tools

Fine-mapping functions are now implemented via
[`echofinemap`](https://github.com/RajLabMSSM/echofinemap):

- `echolocatoR` will automatically check whether you have the necessary
  columns to run each tool you selected in
  `echolocatoR::finemap_loci(finemap_methods=...)`. It will remove any
  tools that for which there are missing necessary columns, and produces
  a message letting you know which columns are missing.
- Note that some columns (e.g. `MAF`,`N`,`t-stat`) will be automatically
  inferred if missing.  
- For easy reference, we list the necessary columns here as well.  
  See `?echodata::construct_colmap()` for descriptions of these
  columns.  
  All methods require the columns: `SNP`,`CHR`,`POS`,`Effect`,`StdErr`

  

``` r
fm_methods <- echofinemap::required_cols(add_versions = FALSE, 
                                         embed_links = TRUE,
                                         verbose = FALSE)
knitr::kable(x = fm_methods)
```

| method | required | suggested | source | citation |
|:---|:---|:---|:---|:---|
| ABF | SNP, CHR…. |  | [source](https://github.com/chr1swallace/coloc) | [cite](https://doi.org/10.1086%2F519024) |
| COJO_conditional | SNP, CHR…. | Freq, P, N | [source](https://github.com/jianyangqt/gcta) | [cite](https://doi.org/10.1038/ng.2213) |
| COJO_joint | SNP, CHR…. | Freq, P, N | [source](https://github.com/jianyangqt/gcta) | [cite](https://doi.org/10.1038/ng.2213) |
| COJO_stepwise | SNP, CHR…. | Freq, P, N | [source](https://github.com/jianyangqt/gcta) | [cite](https://doi.org/10.1038/ng.2213) |
| FINEMAP | SNP, CHR…. | A1, A2, …. | [source](http://www.christianbenner.com/) | [cite](https://doi.org/10.1093%2Fbioinformatics%2Fbtw018) |
| PAINTOR | SNP, CHR…. | MAF | [source](https://github.com/gkichaev/PAINTOR_V3.0) | [cite](https://doi.org/10.1093/bioinformatics/btw615) |
| POLYFUN_FINEMAP | SNP, CHR…. | MAF, N | [source](https://github.com/omerwe/polyfun) | [cite](https://doi.org/10.1038/s41588-022-01036-9) |
| POLYFUN_SUSIE | SNP, CHR…. | MAF, N | [source](https://github.com/omerwe/polyfun) | [cite](https://doi.org/10.1038/s41588-022-01036-9) |
| SUSIE | SNP, CHR…. | N | [source](https://github.com/stephenslab/susieR) | [cite](https://doi.org/10.1371/journal.pgen.1010299) |

## Datasets

Datasets are now stored/retrieved via the following **echoverse**
subpackages:  
- [`echodata`](https://github.com/RajLabMSSM/echodata): Pre-computed
fine-mapping results. Also handles the semi-automated standardization of
summary statistics.  
- [`echoannot`](https://github.com/RajLabMSSM/echoannot): Annotates
GWAS/QTL summary statistics using epigenomics, pre-compiled annotation
matrices, and machine learning model predictions of variant-specific
functional impacts.  
- [`catalogueR`](https://github.com/RajLabMSSM/catalogueR): Large
compendium of fully standardized e/s/t-QTL summary statistics.

For more detailed information about each dataset, use `?`:

``` r
### Examples ###

library(echoannot)   
?NOTT_2019.interactome # epigenomic annotations
library(echodata) 
?BST1 # fine-mapping results 
```

### [**`MungeSumstats`**](https://github.com/neurogenomics/MungeSumstats):

- You can search, import, and standardize any GWAS in the [*Open
  GWAS*](https://gwas.mrcieu.ac.uk/) database via
  [`MungeSumstats`](https://github.com/neurogenomics/MungeSumstats),
  specifically the functions `find_sumstats` and `import_sumstats`.

### [`catalogueR`](https://github.com/RajLabMSSM/catalogueR): QTLs

#### [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/): `catalogueR::eQTL_Catalogue.query()`

- API access to full summary statistics from many standardized e/s/t-QTL
  datasets.  
- Data access and colocalization tests facilitated through the
  [`catalogueR`](https://github.com/RajLabMSSM/catalogueR) R package.

### [`echodata`](https://github.com/RajLabMSSM/catalogueR): fine-mapping results

#### [***echolocatoR Fine-mapping Portal***](https://rajlab.shinyapps.io/Fine_Mapping_Shiny): pre-computed fine-mapping results

- You can visit the *echolocatoR Fine-mapping Portal* to interactively
  visualize and download pre-computed fine-mapping results across a
  variety of phenotypes.
- This data can be searched and imported programmatically using
  [`echodata::portal_query()`](https://rdrr.io/pkg/echodata/man/portal_query.html).

### [`echoannot`](https://github.com/RajLabMSSM/echoannot): Epigenomic & genome-wide annotations

#### [Nott et al. (2019)](https://science.sciencemag.org/content/366/6469/1134.abstract): `echoannot::NOTT2019_*()`

- Data from this publication contains results from cell type-specific
  (neurons, oligodendrocytes, astrocytes, microglia, & peripheral
  myeloid cells) epigenomic assays (H3K27ac, ATAC, H3K4me3) from *ex
  vivo* pediatric human brain tissue.

#### [Corces et al.2020](https://doi.org/10.1038/s41588-020-00721-x): `echoannot::CORCES2020_*()`

- Data from this publication contains results from single-cell and bulk
  chromatin accessibility assays (\[sc\]ATAC-seq) and chromatin
  interactions ( [`FitHiChIP`](https://ay-lab.github.io/FitHiChIP/))
  from *postmortem* adult human brain tissue.

#### [XGR](http://xgr.r-forge.r-project.org): `echoannot::XGR_download_and_standardize()`

- API access to a diverse library of cell type/line-specific epigenomic
  (e.g. **ENCODE**) and other genome-wide annotations.

#### [Roadmap](http://www.roadmapepigenomics.org): `echoannot::ROADMAP_query()`

- API access to cell type-specific epigenomic data.

#### [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html): `echoannot::annotate_snps()`

- API access to various genome-wide SNP annotations (e.g. missense,
  nonsynonmous, intronic, enhancer).

#### [HaploR](https://cran.r-project.org/web/packages/haploR/vignettes/haplor-vignette.html): `echoannot::annotate_snps()`

- API access to known per-SNP QTL and epigenomic data hits.

## Enrichment tools

Annotation enrichment functions are now implemented via
[`echoannot`](https://github.com/RajLabMSSM/echoannot):

### Implemented

#### [XGR](http://xgr.r-forge.r-project.org): `echoannot::XGR_enrichment()`

- Binomial enrichment tests between customisable foreground and
  background SNPs.

#### [motifbreakR](https://github.com/Simon-Coetzee/motifBreakR): `echoannot::MOTIFBREAKR()`

- Identification of transcript factor binding motifs (TFBM) and
  prediction of SNP disruption to said motifs.
- Includes a comprehensive list of TFBM databases via
  [MotifDB](https://bioconductor.org/packages/release/bioc/html/MotifDb.html)
  (9,900+ annotated position frequency matrices from 14 public sources,
  for multiple organisms).

#### [regioneR](http://bioconductor.org/packages/release/bioc/html/regioneR.md): `echoannot::test_enrichment()`

- Iterative pairwise permutation testing of overlap between all
  combinations of two
  [`GRangesList`](https://biodatascience.github.io/compbio/bioc/GRL.html)
  objects.

### Under construction

#### [GARFIELD](https://www.bioconductor.org/packages/release/bioc/html/garfield.html)

- Genomic enrichment with LD-informed heuristics.

#### [GoShifter](https://github.com/immunogenomics/goshifter)

- LD-informed iterative enrichment analysis.

#### [S-LDSC](https://www.nature.com/articles/ng.3954)

- Genome-wide stratified LD score regression.
- Inlccles 187-annotation baseline model from [Gazal et al.
  2018](https://www.nature.com/articles/s41588-018-0231-8).  
- You can alternatively supply a custom annotations matrix.

## LD reference panels

LD reference panels are now queried/processed by
[`echoLD`](https://github.com/RajLabMSSM/echoLD), specifically the
function `get_LD()`:

### [UK Biobank](https://www.ukbiobank.ac.uk)

### [1000 Genomes Phase 1](https://www.internationalgenome.org)

### [1000 Genomes Phase 3](https://www.internationalgenome.org)

### Custom LD panel:

- From user-supplied VCFs

### Custom LD panel

- From user-supplied precomputed LD matrices

## Plotting

Plotting functions are now implemented via:  
- [`echoplot`](https://github.com/RajLabMSSM/echoplot): Multi-track
locus plots with GWAS, fine-mapping results, and functional annotations
(`plot_locus()`). Can also plot multi-GWAS/QTL and multi-ancestry
results (`plot_locus_multi()`).  
- [`echoannot`](https://github.com/RajLabMSSM/echoannot): Study-level
summary plots showing aggregted info across many loci at once
([`super_summary_plot()`](https://rdrr.io/pkg/echoannot/man/super_summary_plot.html)).  
- [`echoLD`](https://github.com/RajLabMSSM/echoLD): Plot an LD matrix
using one of several differnt plotting methods (`plot_LD()`).

## Tabix queries

All queries of [`tabix`](http://www.htslib.org/doc/tabix.md)-indexed
files (for rapid data subset extraction) are implemented via
[`echotabix`](https://github.com/RajLabMSSM/echotabix).

- [`echotabix::convert_and_query()`](https://rdrr.io/pkg/echotabix/man/convert_and_query.html)
  detects whether the GWAS summary statistics file you provided is
  already `tabix`-indexed, and it not, automatically performs all steps
  necessary to convert it (sorting, `bgzip`-compression, indexing)
  across a wide variety of scenarios.  
- [`echotabix::query()`](https://rdrr.io/pkg/echotabix/man/query.html)
  contains many different methods for making tabix queries
  (e.g. `Rtracklayer`,`echoconda`,`VariantAnnotation`,`seqminer`), each
  of which fail in certain circumstances. To avoid this, `query()`
  automatically selects the method that will work for the particular
  file being queried and your machine’s particular versions of
  R/Bioconductor/OS, taking the guesswork and troubleshooting out of
  `tabix` queries.

## Downloads

Single- and multi-threaded downloads are now implemented via
[`downloadR`](https://github.com/RajLabMSSM/downloadR).

- Multi-threaded downloading is performed using
  [`axel`](https://github.com/axel-download-accelerator/axel), and is
  particularly useful for speeding up downloads of large files.
- `axel` is installed via the official *echoverse*
  [conda](https://docs.conda.io/en/latest/) environment: “echoR_mini”.
  This environment is automatically created by the function
  [`echoconda::yaml_to_env()`](https://rdrr.io/pkg/echoconda/man/yaml_to_env.html)
  when needed.

------------------------------------------------------------------------

[Brian M. Schilder, Bioinformatician
II](https://bschilder.github.io/BMSchilder/)  
[Raj Lab](https://rajlab.org)  
[Department of Neuroscience, Icahn School of Medicine at Mount
Sinai](https://icahn.mssm.edu/about/departments-offices/neuroscience)

------------------------------------------------------------------------

# Session info

``` r
utils::sessionInfo()
```

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Tahoe 26.3.1
    ##
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    ##
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ##
    ## time zone: America/New_York
    ## tzcode source: internal
    ##
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base
    ##
    ## loaded via a namespace (and not attached):
    ##   [1] splines_4.5.1               aws.s3_0.3.22
    ##   [3] BiocIO_1.20.0               bitops_1.0-9
    ##   [5] filelock_1.0.3              cellranger_1.1.0
    ##   [7] tibble_3.3.1                R.oo_1.27.1
    ##   [9] basilisk.utils_1.22.0       graph_1.88.1
    ##  [11] rpart_4.1.24                XML_3.99-0.22
    ##  [13] lifecycle_1.0.5             httr2_1.2.2
    ##  [15] mixsqp_0.3-54               rprojroot_2.1.1
    ##  [17] OrganismDbi_1.52.0          ensembldb_2.34.0
    ##  [19] lattice_0.22-9              MASS_7.3-65
    ##  [21] backports_1.5.0             magrittr_2.0.4
    ##  [23] Hmisc_5.2-5                 openxlsx_4.2.8.1
    ##  [25] rmarkdown_2.30              yaml_2.3.12
    ##  [27] dlstats_0.1.7               otel_0.2.0
    ##  [29] zip_2.3.3                   ggbio_1.58.0
    ##  [31] reticulate_1.45.0           gld_2.6.8
    ##  [33] DBI_1.3.0                   RColorBrewer_1.1-3
    ##  [35] abind_1.4-8                 expm_1.0-0
    ##  [37] rvcheck_0.2.1               GenomicRanges_1.62.1
    ##  [39] purrr_1.2.1                 R.utils_2.13.0
    ##  [41] AnnotationFilter_1.34.0     biovizBase_1.58.0
    ##  [43] BiocGenerics_0.56.0         RCurl_1.98-1.17
    ##  [45] nnet_7.3-20                 yulab.utils_0.2.4
    ##  [47] VariantAnnotation_1.56.0    rappdirs_0.3.4
    ##  [49] rworkflows_1.0.8            IRanges_2.44.0
    ##  [51] S4Vectors_0.48.0            echoLD_0.99.12
    ##  [53] echofinemap_1.0.0           irlba_2.3.7
    ##  [55] gitcreds_0.1.2              echodata_1.0.0
    ##  [57] piggyback_0.1.5             codetools_0.2-20
    ##  [59] DelayedArray_0.36.0         DT_0.34.0
    ##  [61] xml2_1.5.2                  tidyselect_1.2.1
    ##  [63] UCSC.utils_1.6.1            farver_2.1.2
    ##  [65] viridis_0.6.5               matrixStats_1.5.0
    ##  [67] stats4_4.5.1                base64enc_0.1-6
    ##  [69] Seqinfo_1.0.0               echotabix_1.0.0
    ##  [71] GenomicAlignments_1.46.0    jsonlite_2.0.0
    ##  [73] e1071_1.7-17                Formula_1.2-5
    ##  [75] survival_3.8-6              tools_4.5.1
    ##  [77] DescTools_0.99.60           Rcpp_1.1.1
    ##  [79] glue_1.8.0                  gridExtra_2.3
    ##  [81] SparseArray_1.10.9          xfun_0.56
    ##  [83] here_1.0.2                  MatrixGenerics_1.22.0
    ##  [85] GenomeInfoDb_1.46.2         dplyr_1.2.0
    ##  [87] withr_3.0.2                 BiocManager_1.30.27
    ##  [89] fastmap_1.2.0               basilisk_1.22.0
    ##  [91] boot_1.3-32                 digest_0.6.39
    ##  [93] R6_2.6.1                    colorspace_2.1-2
    ##  [95] dichromat_2.0-0.1           RSQLite_2.4.6
    ##  [97] cigarillo_1.0.0             R.methodsS3_1.8.2
    ##  [99] tidyr_1.3.2                 generics_0.1.4
    ## [101] renv_1.1.8                  data.table_1.18.2.1
    ## [103] rtracklayer_1.70.1          class_7.3-23
    ## [105] httr_1.4.8                  htmlwidgets_1.6.4
    ## [107] S4Arrays_1.10.1             pkgconfig_2.0.3
    ## [109] gtable_0.3.6                Exact_3.3
    ## [111] blob_1.3.0                  S7_0.2.1
    ## [113] XVector_0.50.0              echoconda_1.0.0
    ## [115] htmltools_0.5.9             susieR_0.14.2
    ## [117] RBGL_1.86.0                 ProtGenerics_1.42.0
    ## [119] scales_1.4.0                Biobase_2.70.0
    ## [121] lmom_3.2                    png_0.1-8
    ## [123] knitr_1.51                  rstudioapi_0.18.0
    ## [125] reshape2_1.4.5              tzdb_0.5.0
    ## [127] rjson_0.2.23                checkmate_2.3.4
    ## [129] badger_0.2.5                curl_7.0.0
    ## [131] proxy_0.4-29                cachem_1.1.0
    ## [133] stringr_1.6.0               rootSolve_1.8.2.4
    ## [135] parallel_4.5.1              foreign_0.8-91
    ## [137] AnnotationDbi_1.72.0        restfulr_0.0.16
    ## [139] desc_1.4.3                  pillar_1.11.1
    ## [141] grid_4.5.1                  reshape_0.8.10
    ## [143] vctrs_0.7.1                 cluster_2.1.8.2
    ## [145] htmlTable_2.4.3             evaluate_1.0.5
    ## [147] readr_2.2.0                 GenomicFeatures_1.62.0
    ## [149] mvtnorm_1.3-3               cli_3.6.5
    ## [151] compiler_4.5.1              echogithub_0.99.5
    ## [153] Rsamtools_2.26.0            rlang_1.1.7
    ## [155] crayon_1.5.3                aws.signature_0.6.0
    ## [157] forcats_1.0.1               plyr_1.8.9
    ## [159] fs_1.6.7                    stringi_1.8.7
    ## [161] coloc_5.2.3                 echoannot_1.0.1
    ## [163] viridisLite_0.4.3           BiocParallel_1.44.0
    ## [165] Biostrings_2.78.0           lazyeval_0.2.2
    ## [167] gh_1.5.0                    Matrix_1.7-4
    ## [169] downloadR_1.0.0             dir.expiry_1.18.0
    ## [171] BSgenome_1.78.0             patchwork_1.3.2
    ## [173] hms_1.1.4                   bit64_4.6.0-1
    ## [175] ggplot2_4.0.2               KEGGREST_1.50.0
    ## [177] haven_2.5.5                 SummarizedExperiment_1.40.0
    ## [179] memoise_2.0.1               snpStats_1.60.0
    ## [181] bit_4.6.0                   readxl_1.4.5

  
