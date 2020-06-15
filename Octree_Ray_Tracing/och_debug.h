#pragma once

#include "och_float.h"
#include "och_string_util.h"

#define OCH_FLT_OUT_SET(f, pad_name, pad_float) '\n' << och::pad(#f, pad_name) << ": " << och::pad(std::to_string(f), pad_float) << " | " << och::fbout(f)

#define OCH_FLT_OUT(f) OCH_FLT_OUT_SET(f, 10, 12)

#define OCH_DEBUG 1

#ifdef OCH_DEBUG

#define OCH_IF_DEBUG(arg) arg

#else

#define OCH_IF_DEBUG(arg)

#endif // OCH_DEBUG
