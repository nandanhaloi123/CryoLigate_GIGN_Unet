#!/bin/bash
chimerax_path=/nethome/emankov/software/packages/chimerax/1.9/usr
export PATH=$chimerax_path/lib/ucsf-chimerax/bin:$PATH
export LD_LIBRARY_PATH=$chimerax_path/lib/ucsf-chimerax/lib:$chimerax_path/lib-extra:$LD_LIBRARY_PATH
export MANPATH=$chimerax_path/share/man
alias chimerax=ChimeraX