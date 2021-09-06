#pragma once

#define RETURN_WITH_COND_TRUE(expr) \
    if (expr)                       \
        return 1;