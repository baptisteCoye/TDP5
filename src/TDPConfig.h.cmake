#ifndef CONFIG_H
#define CONFIG_H

#include <string.h>

char testDataPath [] = "@CMAKE_SOURCE_DIR@/tests/data/";

#cmakedefine VERBOSE @PROJET_VERBOSE@

#cmakedefine SAVE_RESULTS

#endif