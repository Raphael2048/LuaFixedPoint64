#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <inttypes.h>
#include <stdlib.h>
#include <ctype.h>

#include "lua.h"
#include "lauxlib.h"
#include "sinlut.h"


typedef int64_t fp64_t;
static const fp64_t fp64_one = 0x0000000100000000;
static const fp64_t fp64_pi  = 0x3243F6A88;
static const fp64_t fp64_e = 0x2B7E15163;
static const fp64_t fp64_piover2 = 0x1921FB544;
static const fp64_t fp64_pitimes2 = 0x6487ED511;
static const fp64_t fp64_rad2deg = 0x394BB834D1;
static const fp64_t fp64_deg2rad = 0x477D1A8;
static const fp64_t fp64_epsilon = 0x100;
static const fp64_t fp64_zero = 0;
static const fp64_t fp64_piover4 = 0xC90FDAA2;
static const fp64_t fp64_3piover4 = 0x25b2f8fe6;
static const fp64_t fp64_max = 0x7FFFFFFFFFFFFFFF;
static const fp64_t fp64_min = 0x8000000000000000;


#ifdef USELIGHTUSERDATA
#if INTPTR_MAX != INT64_MAX
#undef USELIGHTUSERDATA
#endif
#endif

static fp64_t fpmul(fp64_t a, fp64_t b)
{
    int64_t A = a >> 32, C = b >> 32;
    uint64_t B = a & 0xFFFFFFFF, D = b & 0xFFFFFFFF;
    int64_t AC = A * C;
    int64_t AD_CB = A * D + C * B;
    uint64_t BD = B * D;

    int64_t product_hi = AC + (AD_CB >> 32);

    uint64_t ad_cb_temp = AD_CB << 32;
    uint64_t product_lo = BD + ad_cb_temp;

    if (product_lo < BD)
        product_hi++;

    uint64_t product_lo_tmp = product_lo;
    product_lo -= 0x80000000;
    product_lo -= (uint64_t)product_hi >> 63;
    if (product_lo > product_lo_tmp)
        product_hi--;
    fp64_t result = (product_hi << 32) | (product_lo >> 32);
    result += 1;
    return result;
}

#ifdef __GNUC__
#define clz(x) (__builtin_clzll(x) - (8 * sizeof(long long) - 64))
#else
static uint8_t clz(uint64_t x)
{
    uint8_t result = 0;
    if (x == 0) return 32;
    while (!(x & 0xF000000000000000)) { result += 4; x <<= 4; }
    while (!(x & 0x8000000000000000)) { result += 1; x <<= 1; }
    return result;
}
#endif

static fp64_t fpdiv(fp64_t a, fp64_t b)
{
    if (b == 0) return 0;
    uint64_t remainder = (a >= 0) ? a : -a;
    uint64_t divider = (b >= 0) ? b : -b;
    uint64_t quotient = 0;
    int bit_pos = 33;

    while (!(divider & 0xF) && bit_pos >= 4)
    {
        divider >>= 4;
        bit_pos -= 4;
    }

    while (remainder && bit_pos >= 0)
    {
        int shift = clz(remainder);
        if (shift > bit_pos) shift = bit_pos;
        remainder <<= shift;
        bit_pos -= shift;

        uint64_t div = remainder / divider;
        remainder = remainder % divider;
        quotient += div << bit_pos;

        remainder <<= 1;
        bit_pos--;
    }
    quotient++;
    fp64_t result = quotient >> 1;
    if ((a ^ b) & 0x8000000000000000)
    {
        result = -result;
    }
    return result;
}

static fp64_t fpsqrt(fp64_t value)
{
    uint8_t neg = (value < 0);
    uint64_t num = (neg ? -value : value);
    uint64_t result = 0;
    uint64_t bit;
    uint8_t n;

    //most values <= 65535
    if (num & 0xFFFF000000000000)
        bit = (uint64_t)1 << 62;
    else
        bit = (uint64_t)1 << 46;
    while (bit > num) bit >>= 2;

    for (n = 0; n < 2; n++)
    {
        while (bit)
        {
            if (num >= result + bit)
            {
                num -= result + bit;
                result = (result >> 1) + bit;
            }
            else
            {
                result = result >> 1;
            }
            bit >>= 2;
        }

        if (n == 0)
        {
            //frac last 16bits
            if (num > 0xFFFFFFFF)
            {
                num -= result;
                num = (num << 32) - 0x80000000;
                result = (result << 32) + 0x80000000;
            }
            else
            {
                num <<= 32;
                result <<= 32;
            }
            bit = 1 << 30;
        }
    }

    //rounding
    if (num > result)
    {
        result++;
    }

    return (neg ? -(fp64_t)result : (fp64_t)result);
}

static fp64_t fpsin(fp64_t value)
{
    fp64_t temp = value % (fp64_pitimes2);
    if (temp < 0) temp += fp64_pitimes2;
    bool  sign = true;
    if (temp >= fp64_pi) {
        temp -= fp64_pi;
        sign = false;   
    }
    if (temp >= fp64_piover2)
        temp = fp64_pi - temp;

    uint32_t index = temp >> 15;
    
    fp64_t result = fp64_one;
    if (index <= sin_lut_count) {
        fp64_t frac = (temp & 0x7FFF) << 17;
        fp64_t a = sin_lut[index];
        fp64_t b = index == sin_lut_count ? fp64_one : sin_lut[index + 1];
        result = fpmul(a, fp64_one - frac) + fpmul(b, frac);
    }
    return sign ? result : -result;
}

static fp64_t fpcos(fp64_t value)
{
    return fpsin(value + fp64_piover2);
}

static fp64_t fptan(fp64_t value)
{
    return fpdiv(fpsin(value), fpcos(value));
}

static fp64_t fpatan2(fp64_t y, fp64_t x)
{
    fp64_t abs_inY, mask, angle, r, r_3;
    mask = (y >> (sizeof(fp64_t) * 8 - 1));
    abs_inY = (y + mask) ^ mask;

    if (x >= 0) {
        r = fpdiv(x - abs_inY, x + abs_inY);
        r_3 = fpmul(fpmul(r, r), r);
        angle = fpmul(0x32400000, r_3) - fpmul(0xFB500000, r) + fp64_piover4;
    }
    else {
        r = fpdiv(x + abs_inY, abs_inY - x);
        r_3 = fpmul(fpmul(r, r), r);
        angle = fpmul(0x32400000, r_3) - fpmul(0xFB500000, r) + fp64_3piover4;
    }
    return y > 0 ? angle : -angle;
}

static fp64_t fpatan(fp64_t v)
{
    return fpatan2(v, fp64_one);
}

static fp64_t fpasin(fp64_t x)
{
    if (x > fp64_one || x < -fp64_one)
    {
        return fp64_zero;
    }
    fp64_t out;
    out = (fp64_one - fpmul(x, x));
    out = fpdiv(x, fpsqrt(out));
    out = fpatan(out);
    return out;
}

static fp64_t fpacos(fp64_t x)
{
    return (fp64_piover2 - fpasin(x));
}

static inline fp64_t from_number(lua_Number l)
{
    lua_Number temp = l * fp64_one;
    temp += (temp >= 0) ? 0.5f : -0.5f;
    return (fp64_t)(temp);
}
static inline fp64_t from_integer(lua_Integer d)
{ 
    return d * fp64_one;
}

static fp64_t fpexp(fp64_t x)
{
    if (x == 0) return fp64_one;
    if (x == fp64_one) return fp64_e;
    if (x >= 92288378626) return fp64_max;
    if (x <= -95265423098) return fp64_zero;

    bool neg = x < 0;
    if (neg) x = -x;

    fp64_t result = x + fp64_one;
    fp64_t term = x;

    for (int i = 2; i < 30; i++)
    {
        term = fpmul(term, fpdiv(x, from_integer(i)));
        result += term;

        if ((term < 500000) && ((i > 15) || (term < 20000)))
            break;
    }
    if (neg) result = fpdiv(fp64_one, result);
    return result;
}

static fp64_t fplog(fp64_t x)
{
    fp64_t guess = from_integer(2);
    fp64_t delta;
    int scaling = 0;
    int count = 0;

    if (x <= 0)
        return fp64_min;

    // Bring the value to the most accurate range (1 < x < 100)
    const fp64_t e_to_fourth = 234497268814;
    while (x > from_integer(2))
    {
        x = fpdiv(x, e_to_fourth);
        scaling += 4;
    }

    while (x < fp64_one)
    {
        x = fpmul(x, e_to_fourth);
        scaling -= 4;
    }

    do
    {
        // Solving e(x) = y using Newton's method
        // f(x) = e(x) - y
        // f'(x) = e(x)
        fp64_t e = fpexp(guess);
        delta = fpdiv(x - e, e);

        // It's unlikely that logarithm is very large, so avoid overshooting.
        if (delta > from_integer(3))
            delta = from_integer(3);

        guess += delta;
    } while ((count++ < 10)
        && ((delta > 1) || (delta < -1)));

    return guess + from_integer(scaling);
}

static inline fp64_t fprs(fp64_t x)
{
    fp64_t y = (x >> 1) + (x & 1);
    return y;
}

static fp64_t fplog2_inner(fp64_t x)
{
    fp64_t result = 0;

    while (x >= from_integer(2))
    {
        result++;
        x = fprs(x);
    }

    if (x == 0) return (result << 32);

    uint_fast8_t i;
    for (i = 32; i > 0; i--)
    {
        x = fpmul(x, x);
        result <<= 1;
        if (x >= from_integer(2))
        {
            result |= 1;
            x = fprs(x);
        }
    }
    x = fpmul(x, x);
    if (x >= from_integer(2)) result++;

    return result;
}

fp64_t fplog2(fp64_t x)
{
    if (x <= 0) return fp64_min;

    // If the input is less than one, the result is -log2(1.0 / in)
    if (x < fp64_one)
    {
        fp64_t inverse = fpdiv(fp64_one, x);
        return -fplog2_inner(inverse);
    }

    // If input >= 1, just proceed as normal.
    // Note that x == fix16_one is a special case, where the answer is 0.
    return fplog2_inner(x);
}

static fp64_t from_str(const char* buf)
{
    while (isspace(*buf))
        buf++;
    bool negative = (*buf == '-');
    if (*buf == '+' || *buf == '-')
        buf++;

    uint64_t intpart = 0;
    int count = 0;
    while (isdigit(*buf))
    {
        intpart *= 10;
        intpart += *buf++ - '0';
        count++;
    }
    fp64_t value = intpart << 32;

    if (*buf == '.' || *buf == ',')
    {
        buf++;

        int64_t fracpart = 0;
        int64_t scale = 1;
        while (isdigit(*buf) && scale < 10000000000)
        {
            scale *= 10;
            fracpart *= 10;
            fracpart += *buf++ - '0';
        }

        value = value + fpdiv(fracpart, scale);
    }
    return negative ? -value : value;
}

static char* itoa_loop(char* buf, uint64_t scale, uint64_t value, bool skip)
{
    while (scale)
    {
        unsigned digit = (value / scale);
        if (!skip || digit || scale == 1)
        {
            skip = false;
            *buf++ = '0' + digit;
            value %= scale;
        }
        scale /= 10;
    }
    return buf;
}

static void to_str(fp64_t value, char* buf)
{
    uint64_t uvalue = (value > 0) ? value : -value;
    if (value < 0)
        *buf++ = '-';
    uint64_t intpart = uvalue >> 32;
    uint64_t fracpart = uvalue & 0xFFFFFFFF;
    uint64_t scale = 1000000000;
    fracpart = fpmul(fracpart, scale);
    if (fracpart >= scale)
    {
        intpart++;
        fracpart -= scale;
    }
    buf = itoa_loop(buf, 1000000000, intpart, true);
    if (fracpart != 0)
    {
        *buf++ = '.';
        buf = itoa_loop(buf, scale / 10, fracpart, false);
        buf--;
        while(*buf == '0') buf--;
        buf++;
    }
    *buf = '\0';
}

static inline fp64_t fpceil(fp64_t n) { return (n & 0xFFFFFFFF00000000) + (n & 0x00000000FFFFFFFF ? fp64_one : 0); }
static inline fp64_t fpfloor(fp64_t n) { return (n & 0xFFFFFFFF00000000); }
static inline fp64_t fpabs(fp64_t x) { return (x < 0 ? -x : x); }
static inline fp64_t fpmax(fp64_t x, fp64_t y) { return x > y ? x : y; }
static inline fp64_t fpmin(fp64_t x, fp64_t y) { return x < y ? x : y; }
static inline fp64_t fpclamp(fp64_t x, fp64_t y, fp64_t z) { return fpmin(fpmax(x, y), z);}

static inline lua_Number to_number(fp64_t n) { return (lua_Number)n / fp64_one; }

static bool _isfp64(lua_State* L, int pos)
{
    if (lua_getmetatable(L, pos))
    {
        luaL_getmetatable(L, "fp64");
        int equal = lua_rawequal(L, -1, -2);
        lua_pop(L, 2);
        return equal;
    }
    return false;
}

LUALIB_API void pushfp64(lua_State* L, fp64_t n)
{
#ifdef USELIGHTUSERDATA
    lua_pushlightuserdata(L, (void*)n);
#else
    fp64_t* p = (fp64_t*)lua_newuserdata(L, sizeof(fp64_t));
    *p = n;
    luaL_getmetatable(L, "fp64");
    lua_setmetatable(L, -2);
#endif
}

static inline fp64_t _parse_str(lua_State* L, int pos)
{
    const char* str = lua_tostring(L, pos);
    return from_str(str);
}

static int64_t tofp64(lua_State* L, int pos)
{
    fp64_t n = 0;
    int type = lua_type(L, pos);

    switch(type) 
    {
        case LUA_TNUMBER:
#if !defined(LUA_VERSION_NUM) || LUA_VERSION_NUM < 503
            if (lua_tonumber(L, pos) == lua_tointeger(L, pos))
#else
            if (lua_isinteger(L, pos))
#endif
            {
                lua_Integer l = lua_tointeger(L, pos);
                if (l == 0)
                {
                    n = fp64_zero;
                }
                else
                {
                    n = from_integer(lua_tointeger(L, pos));
                }
            }
            else
            {
                n = from_number(lua_tonumber(L, pos));
            }
            break;
        case LUA_TSTRING:
            n = _parse_str(L, n);
            break;
#ifdef USELIGHTUSERDATA
        case LUA_TLIGHTUSERDATA:
            n = (int64_t)lua_touserdata(L, pos);
#else
        case LUA_TUSERDATA:
            if (_isfp64(L, pos))
            {
                n = *(fp64_t*)lua_touserdata(L, pos);            
            }
#endif 
            break;
        default:
            n = fp64_zero;
            const char *msg = lua_pushfstring(L, "%s expected, got %s", "fp64", luaL_typename(L, pos));
            return luaL_argerror(L, pos, msg);        
            break;
    }
    return n;
}

static int _fp64add(lua_State* L)
{
    fp64_t lhs = tofp64(L, 1);    
    fp64_t rhs = tofp64(L, 2);
    pushfp64(L, lhs + rhs);
    return 1;
}

static int _fp64sub(lua_State* L)
{
    fp64_t lhs = tofp64(L, 1);    
    fp64_t rhs = tofp64(L, 2);
    pushfp64(L, lhs - rhs);
    return 1;
}

static int _fp64mul(lua_State* L)
{
    fp64_t lhs = tofp64(L, 1);    
    fp64_t rhs = tofp64(L, 2);
    pushfp64(L, fpmul(lhs, rhs));
    return 1;    
}

static int _fp64div(lua_State* L)
{
    fp64_t lhs = tofp64(L, 1);    
    fp64_t rhs = tofp64(L, 2);

    if (rhs == 0) 
    {
        return luaL_error(L, "div by zero");
    }

    pushfp64(L, fpdiv(lhs , rhs));
    return 1;
}

static int _fp64mod(lua_State* L)
{
    fp64_t lhs = tofp64(L, 1);    
    fp64_t rhs = tofp64(L, 2);
    if (rhs == 0) 
    {
        return luaL_error(L, "mod by zero");
    }
    pushfp64(L, lhs % rhs);
    return 1;
}

static int _fp64unm(lua_State* L)
{
    fp64_t lhs = *(fp64_t*)lua_touserdata(L, 1);        
    pushfp64(L, -lhs);
    return 1;
}

static int _fp64max(lua_State * L)
{
    fp64_t lhs = tofp64(L, 1);    
    fp64_t rhs = tofp64(L, 2);
    pushfp64(L, fpmax(lhs, rhs));
    return 1;
}

static int _fp64min(lua_State * L)
{
    fp64_t lhs = tofp64(L, 1);    
    fp64_t rhs = tofp64(L, 2);
    pushfp64(L, fpmin(lhs, rhs));
    return 1;
}

static int _fp64clamp(lua_State * L)
{
    fp64_t x = tofp64(L, 1);    
    fp64_t lo = tofp64(L, 2);
    fp64_t hi = tofp64(L, 3);
    pushfp64(L, fpclamp(x, lo, hi));
    return 1;
}

static int _fp64eq(lua_State* L)
{
    fp64_t lhs = *(fp64_t*)lua_touserdata(L, 1);     
    fp64_t rhs = *(fp64_t*)lua_touserdata(L, 2);     
    lua_pushboolean(L, lhs == rhs);    
    return 1;
}

static int _fp64lt(lua_State* L)
{
    fp64_t lhs = tofp64(L, 1); 
    fp64_t rhs = tofp64(L, 2);
    lua_pushboolean(L, lhs < rhs);
    return 1;
}

static int _fp64le(lua_State* L)
{
    fp64_t lhs = tofp64(L, 1); 
    fp64_t rhs = tofp64(L, 2);
    lua_pushboolean(L, lhs <= rhs);
    return 1;
}

static int _fp64tostring(lua_State* L)
{
    fp64_t n = tofp64(L, 1);    
    char temp[32];
    to_str(n, temp);
    lua_pushstring(L, temp);
    return 1;
}

static int _fp64tohex(lua_State* L)
{
    fp64_t n = tofp64(L, 1);    
    char temp[32];
    sprintf(temp, "%llx", n);
    lua_pushstring(L, temp);
    return 1;
}


static int _fp64equals(lua_State* L)
{
    int64_t lhs = tofp64(L, 1);
    int64_t rhs = tofp64(L, 2);
    lua_pushboolean(L, lhs == rhs);
    return 1;
}

static int _fp64compare(lua_State *L)
{
    int64_t lhs = tofp64(L, 1);
    int64_t rhs = tofp64(L, 2);
    int res = lhs == rhs ? 0 : (lhs < rhs ? -1 : 1);
    lua_pushinteger(L, res);
    return 1;
}

static int _fp64tonumber(lua_State *L)
{
    fp64_t n = tofp64(L, 1);
    lua_Number l = to_number(n);
    lua_pushnumber(L, l);
    return 1;
}

static int _fp64sqrt(lua_State* L)
{
    fp64_t n = tofp64(L, 1);
    fp64_t l = fpsqrt(n);
    pushfp64(L, l);
    return 1;
}

static int _fp64ceil(lua_State* L)
{
    fp64_t n = tofp64(L, 1);
    fp64_t l = fpceil(n);
    pushfp64(L, l);
    return 1;
}

static int _fp64floor(lua_State* L)
{
    fp64_t n = tofp64(L, 1);
    fp64_t l = fpfloor(n);
    pushfp64(L, l);
    return 1;
}

static int _fp64abs(lua_State* L)
{
    fp64_t n = tofp64(L, 1);
    fp64_t l = fpabs(n);
    pushfp64(L, l);
    return 1;
}


static int _fp64sin(lua_State* L)
{
    fp64_t angle = tofp64(L, 1);
    pushfp64(L, fpsin(angle));
    return 1;
}

static int _fp64cos(lua_State* L)
{
    fp64_t angle = tofp64(L, 1);
    pushfp64(L, fpcos(angle));
    return 1;
}

static int _fp64tan(lua_State* L)
{
    fp64_t angle = tofp64(L, 1);
    pushfp64(L, fptan(angle));
    return 1;
}

static int _fp64asin(lua_State* L)
{
    fp64_t value = tofp64(L, 1);
    pushfp64(L, fpasin(value));
    return 1;
}
static int _fp64acos(lua_State* L)
{
    fp64_t value = tofp64(L, 1);    
    pushfp64(L, fpacos(value));
    return 1;
}

static int _fp64atan(lua_State* L)
{
    fp64_t value = tofp64(L, 1);    
    pushfp64(L, fpatan(value));
    return 1;
}

static int _fp64atan2(lua_State* L)
{
    fp64_t lhs = tofp64(L, 1);    
    fp64_t rhs = tofp64(L, 2);
    pushfp64(L, fpatan2(lhs, rhs));
    return 1;
}

static int _fp64exp(lua_State* L)
{
    fp64_t x = tofp64(L, 1);    
    pushfp64(L, fpexp(x));
    return 1;
}

static int _fp64log(lua_State* L)
{
    fp64_t x = tofp64(L, 1);    
    pushfp64(L, fplog(x));
    return 1;
}

static int _fp64log2(lua_State* L)
{
    fp64_t x = tofp64(L, 1);    
    pushfp64(L, fplog2(x));
    return 1;
}


static int newfp64(lua_State* L)
{
    int64_t n = fp64_zero;
    int type = lua_type(L, 1);

    if (type == LUA_TSTRING)
    {
        n = _parse_str(L, 1);
    }
    else if (type == LUA_TNUMBER)
    {
#if !defined(LUA_VERSION_NUM) || LUA_VERSION_NUM < 503
        if(lua_tonumber(L, 1) == lua_tointeger(L, 1))
#else
        if(lua_isinteger(L, 1))
#endif
        {
            n = from_integer(lua_tointeger(L, 1));
        }
        else 
        {
            n = from_number(lua_tonumber(L, 1));
        }
    }
    pushfp64(L, n);
    return 1;
}

static const struct luaL_Reg lib_fp64_meta [] = {
    {"__add", _fp64add},
    {"__sub", _fp64sub},
    {"__mul", _fp64mul},
    {"__div", _fp64div},
    {"__mod", _fp64mod},
    {"__unm", _fp64unm},
    {"__eq", _fp64eq},
    {"__lt", _fp64lt},
    {"__le", _fp64le},
    {"__tostring", _fp64tostring},
    {"__index", NULL},
    {NULL, NULL}
};

static const struct luaL_Reg lib_fp64 [] = {
    {"hex", _fp64tohex},
    {"tostring", _fp64tostring},
    {"tonumber", _fp64tonumber},
    {"compare", _fp64compare},
    {"min", _fp64min},
    {"max", _fp64max},
    {"clamp", _fp64clamp},
    {"equals", _fp64equals},
    {"sqrt", _fp64sqrt},
    {"ceil", _fp64ceil},
    {"floor", _fp64floor},
    {"abs", _fp64abs},
    {"sin", _fp64sin},
    {"cos", _fp64cos},
    {"tan", _fp64tan},
    {"asin", _fp64asin},
    {"acos", _fp64acos},
    {"atan", _fp64atan},
    {"atan2", _fp64atan2},
    {"exp", _fp64exp},
    {"log", _fp64log},
    {"log2", _fp64log2},
    {"new", newfp64},
    {"one", NULL},
    {"pi", NULL},
    {"epsilon", NULL},
    {"rad2deg", NULL},
    {"deg2rad", NULL},
    {"e", NULL},
    {NULL, NULL}
};

int luaopen_fp64(lua_State* L)
{
    luaL_newmetatable(L, "fp64");
    luaL_setfuncs(L, lib_fp64_meta, 0);
    luaL_newlib(L, lib_fp64);
    pushfp64(L, fp64_pi);
    lua_setfield(L, -2, "pi");
    pushfp64(L, fp64_one);
#ifdef USELIGHTUSERDATA
    //set metatable for light userdata
    lua_pushvalue(L, -3);
    lua_setmetatable(L, -2);
#endif
    lua_setfield(L, -2, "one");
    pushfp64(L, fp64_epsilon);
    lua_setfield(L, -2, "epsilon");
    pushfp64(L, fp64_rad2deg);
    lua_setfield(L, -2, "rad2deg");
    pushfp64(L, fp64_deg2rad);
    lua_setfield(L, -2, "deg2rad");
    pushfp64(L, fp64_e);
    lua_setfield(L, -2, "e");
    lua_pushvalue(L, -1);
    lua_setfield(L, -3, "__index");
    return 1;
}

#ifdef __cplusplus
}
#endif