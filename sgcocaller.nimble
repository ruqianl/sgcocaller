# Package

version       = "0.3.3"
author        = "Ruqian Lyu"
description   = "sgcocaller: calling crossovers from single-gamete DNA sequencing datasets"
license       = "MIT"
srcDir        = "src"
bin           = @["sgcocaller"]


# Dependencies

requires "nim >= 1.0.0","docopt","hts >= 0.3.4"
requires "https://github.com/ruqianl/distributions.git"

# task test, "run the tests":
#   exec "nim c  -d:useSysAssert -d:useGcAssert --lineDir:on --debuginfo -r tests/test_groups"
#   exec "bash tests/functional-tests.sh"
