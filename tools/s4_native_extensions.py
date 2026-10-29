"""Checked mechanical hooks for a separate S4 stripe/rectangular-basis build."""
import difflib
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]

LUA_METHODS = r'''
static int S4L_Simulation_SetRectangularBasis(lua_State *L){
    Simulation *S = S4L_get_simulation(L, 1);
    int mx = luaL_checkinteger(L, 2), my = luaL_checkinteger(L, 3);
    luaL_argcheck(L, S != NULL, 1, "Simulation object required");
    luaL_argcheck(L, mx >= 0 && mx <= 1000, 2, "Invalid x half-width");
    luaL_argcheck(L, my >= 0 && my <= 1000, 3, "Invalid y half-width");
    int count = (2*mx+1)*(2*my+1);
    luaL_argcheck(L, count <= 8191, 2, "Rectangular basis exceeds adapter limit");
    luaL_argcheck(L, S->Lr[2] != 0 || S->Lr[3] != 0, 1, "2D lattice required");
    Simulation_SetNumG(S, count);
    S->options.waferscope_rectangular[0] = mx;
    S->options.waferscope_rectangular[1] = my;
    return 0;
}
static int S4L_Simulation_UseAnalyticStripeRule(lua_State *L){
    Simulation *S = S4L_get_simulation(L, 1);
    luaL_argcheck(L, S != NULL, 1, "Simulation object required");
    luaL_argcheck(L, S->solution == NULL && S->field_cache == NULL, 1,
        "Configure stripe rule before querying fields or basis");
    S->options.waferscope_stripes = lua_gettop(L) < 2 || lua_toboolean(L, 2);
    return 0;
}
'''

RECTANGULAR = '''if(S->options.waferscope_rectangular[0] >= 0){
        int mx = S->options.waferscope_rectangular[0];
        int my = S->options.waferscope_rectangular[1];
        int next = 1;
        sol->G[0] = sol->G[1] = 0;
        for(int x = -mx; x <= mx; ++x){
            for(int y = -my; y <= my; ++y){
                if(x == 0 && y == 0){ continue; }
                sol->G[2*next] = x; sol->G[2*next+1] = y; ++next;
            }
        }
    }else if(0 != S->Lr[2] || 0 != S->Lr[3]){'''


def apply_extensions(source, patch_path):
    source = Path(source)
    changes = {
        'S4/S4.h': [('} Options;', 'int waferscope_stripes;\nint waferscope_rectangular[2];\n} Options;')],
        'S4/S4.cpp': [
            ('S->options.use_experimental_fmm = 0;', 'S->options.use_experimental_fmm = 0;\n'
             'S->options.waferscope_stripes = 0;\nS->options.waferscope_rectangular[0] = -1;\n'
             'S->options.waferscope_rectangular[1] = -1;'),
            ('if(0 != S->Lr[2] || 0 != S->Lr[3]){\n\t\tunsigned int NG',
             RECTANGULAR+'\n\t\tunsigned int NG'),
            ('if(S->options.use_experimental_fmm){',
             'const int ws = S->options.waferscope_stripes ? '
             'WS4StripeMaterial(S,L,n,pB->Epsilon2,pB->Epsilon_inv) : 0;\n'
             'if(ws < 0){ return 1; }\n'
             'if(ws > 0){ /* Exact stripe layer handled by project adapter. */ }\n'
             'else if(S->options.use_experimental_fmm){')],
        'S4/main_lua.c': [
            ('static int S4L_Simulation_SetNumG(lua_State *L){', LUA_METHODS+'\nstatic int S4L_Simulation_SetNumG(lua_State *L){'),
            ('Simulation_SetNumG(S, n);', 'S->options.waferscope_rectangular[0] = -1;\n'
             'S->options.waferscope_rectangular[1] = -1;\nSimulation_SetNumG(S, n);'),
            ('{"SetNumG", S4L_Simulation_SetNumG},', '{"SetNumG", S4L_Simulation_SetNumG},\n'
             '{"SetRectangularBasis", S4L_Simulation_SetRectangularBasis},\n'
             '{"UseAnalyticStripeRule", S4L_Simulation_UseAnalyticStripeRule},')],
    }
    diffs, pending = [], []
    for relative, replacements in changes.items():
        file = source/relative
        original = file.read_text(encoding='utf-8')
        updated = original
        if relative == 'S4/S4.cpp':
            updated = '#include "waferscope_fmm.h"\n'+updated
        for before, after in replacements:
            if updated.count(before) != 1:
                raise ValueError('Unexpected pinned source anchor: '+relative+' / '+before[:60])
            updated = updated.replace(before, after)
        diffs.extend(difflib.unified_diff(original.splitlines(True), updated.splitlines(True),
            fromfile='a/'+relative, tofile='b/'+relative))
        pending.append((file, updated))
    for file, updated in pending:
        file.write_text(updated, encoding='utf-8')
    Path(patch_path).write_text(''.join(diffs), encoding='utf-8')
    header = ROOT/'s4_prototype/native/waferscope_fmm.h'
    shutil.copy2(header, source/'S4/waferscope_fmm.h')
    return header
