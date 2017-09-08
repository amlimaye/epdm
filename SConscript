#bring in environment prepared by sconsutils
Import('env')

#bring in CPPPATH from environment, don't track garden dependencies
cpp = []
flg = []
for p in env['CPPPATH']:
    if p.startswith('/proj') or p.startswith('/gdn'):
        flg.append('-I%s' % p)
    else:
        cpp.append(p)
env.Replace(CPPPATH=cpp)
env.Append(CFLAGS=flg,CXXFLAGS=flg)

#either compile a debug build or a -O3 build
debug = ARGUMENTS.get('debug',0)
if int(debug):
    env.Append(CXXFLAGS='-ggdb -D_LIBCPP_INLINE_VISIBILITY=\"\" -D\"_LIBCPP_EXTERN_TEMPLATE(...)=\" -D__DEBUG')
else:
    env.Append(CXXFLAGS='-O3 -ffast-math')

#bring in the jsoncpp library
env.Append(LIBS='jsoncpp')
env.Program('test.cpp')
