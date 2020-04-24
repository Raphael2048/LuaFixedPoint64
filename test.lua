fp64 = require("fp64")

local f = fp64.new(5.259451158)
f= f / 5
f = -f
f = f * 3
f= f:sin()
print(f + fp64.pi)
print(fp64.tan(6))
print(f:tonumber())
print(f:hex())
local s = fp64.sqrt(fp64.e)
print(s:log())
print(fp64.log2(2 ^ 10))