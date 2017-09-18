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
folding = ARGUMENTS.get('folding',0)
if int(debug):
    if int(folding):
        env.Append(CXXFLAGS='-ggdb -O0 -D_LIBCPP_INLINE_VISIBILITY=\"\" -D\"_LIBCPP_EXTERN_TEMPLATE(...)=\" -D__DEBUG -D__INCLUDE_FOLDING')
    else: 
        env.Append(CXXFLAGS='-ggdb -O0 -D_LIBCPP_INLINE_VISIBILITY=\"\" -D\"_LIBCPP_EXTERN_TEMPLATE(...)=\" -D__DEBUG')
else:
    if int(folding):
        env.Append(CXXFLAGS='-O3 -ffast-math -D__INCLUDE_FOLDING')
    else:
        env.Append(CXXFLAGS='-O3 -ffast-math')

#bring in the jsoncpp library
env.Append(LIBS='jsoncpp')

#construct the executable
src_files = ['species_utilities.cpp','pop_utilities.cpp',
             'rxn_utilities.cpp','ensemble_utilities.cpp',
             'folding.cpp','toplevel.cpp','types.cpp']
src_files_abspaths = map(lambda x: 'src/' + x, src_files)

env.Program('run_epdm',src_files_abspaths)