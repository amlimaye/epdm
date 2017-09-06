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

#some debug-related flags; libc++11 inlines fcns i'd like to call in the debugger
env.Append(CXXFLAGS='-g -D_LIBCPP_INLINE_VISIBILITY=\"\" -D\"_LIBCPP_EXTERN_TEMPLATE(...)=\"')

#bring in the jsoncpp library
env.Append(LIBS='jsoncpp')
env.Program('test.cpp')
