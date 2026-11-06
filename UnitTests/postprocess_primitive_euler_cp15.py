#!/usr/bin/env python3
"""Read-only CFD-data diagnostic: common half-face Cp and axial-force comparison.

Creates a NEW output directory; never runs a solver or modifies archived data.
The primitive 9000-step field is not certified stationary. Native means the
native Float64 field evaluated through primitive reconstruction at step zero.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import numpy as np
from convert_primitive_euler_initial import mesh_coordinates
from run_primitive_euler_driver_e2e import binary
from run_primitive_euler_naca15 import ROOT,CASES,MESH,write
from run_primitive_euler_steady15c import history

PINF,QINF,AREA=1.,.028,1.
REGIONS=['all','LE','central','TE']


def region(x): return 'LE' if x<.1 else ('TE' if x>.9 else 'central')


def surface(path):
    with path.open() as stream: rows=list(csv.DictReader(stream))
    return [dict(marker=r['marker'],node=int(r['node']),**{k:float(r[k]) for k in ['x','y','p','Cp','nx','ny']}) for r in rows]


def geometric_faces():
    lines=MESH.read_text().splitlines();xy=mesh_coordinates(MESH)
    at=next(i for i,x in enumerate(lines) if x.startswith('NELEM='));n=int(lines[at].split('=')[1])
    cells={}
    for line in lines[at+1:at+1+n]:
        values=[int(v) for v in line.split()];assert values[0]==5
        tri=values[1:4];center=np.mean(xy[tri],axis=0)
        for i in range(3):cells.setdefault(tuple(sorted((tri[i],tri[(i+1)%3]))),[]).append(center)
    at=next(i for i,x in enumerate(lines) if x.startswith('MARKER_TAG=') and x.split('=')[1].strip()=='airfoil')
    n=int(lines[at+1].split('=')[1]);faces=[]
    for segment,line in enumerate(lines[at+2:at+2+n]):
        kind,a,b=[int(v) for v in line.split()];assert kind==3 and a!=b
        adjacent=cells[tuple(sorted((a,b)))];assert len(adjacent)==1
        delta=xy[b]-xy[a];normal=np.array([delta[1],-delta[0]])/2
        if np.dot(normal,adjacent[0]-(xy[a]+xy[b])/2)>0:normal=-normal
        for owner,other in [(a,b),(b,a)]:
            q=.75*xy[owner]+.25*xy[other]
            faces.append(dict(segment=segment,node=owner,x=q[0],y=q[1],nx=normal[0],ny=normal[1],
                              side='upper' if q[1]>0 else 'lower',region=region(q[0]),length=float(np.linalg.norm(normal))))
    assert n==200 and len(faces)==400
    assert all(f['y']!=0 for f in faces)
    owners={f['node'] for f in faces}
    assert len(owners)==200 and all(sum(f['node']==i for f in faces)==2 for i in owners)
    np.testing.assert_allclose(np.sum([[f['nx'],f['ny']] for f in faces],axis=0),0,rtol=0,atol=1e-14)
    return faces


def align(geometry,rows):
    assert len(rows)==len(geometry)==400
    used=set();result=[]
    for face in geometry:
        candidates=[i for i,r in enumerate(rows) if r['node']==face['node'] and
                    abs(r['x']-face['x'])<2e-15 and abs(r['y']-face['y'])<2e-15]
        assert len(candidates)==1 and candidates[0] not in used
        i=candidates[0];used.add(i);r=rows[i]
        assert r['marker']=='airfoil' and all(np.isfinite(r[k]) for k in ['x','y','p','Cp','nx','ny'])
        assert r['p']>0
        np.testing.assert_allclose([r['nx'],r['ny']],[face['nx'],face['ny']],rtol=0,atol=2e-15)
        np.testing.assert_allclose(r['Cp'],(r['p']-PINF)/QINF,rtol=2e-14,atol=2e-14)
        result.append(r)
    assert len(used)==400
    assert len({(r['node'],r['x'],r['y']) for r in rows})==400
    assert sum(r['node']==199 for r in rows)==2
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,default=CASES/'results15d/cp_comparison')
    a=p.parse_args();out=a.output_dir.resolve()
    folder=CASES/'results15d/run_final_native_PROJECTED_cfl005'
    sources=[MESH,CASES/'results/native/native_state_000599.dat',CASES/'results15d/native_initial.dat',
             folder/'primitive_state_00000.dat',folder/'primitive_state_09000.dat',folder/'history.csv',
             folder/'surface_primitive_pressure_0.csv',folder/'surface_primitive_pressure_9000.csv',
             CASES/'results/native/history.csv']
    assert all(path.is_file() for path in sources)
    hashes={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
    original=binary(sources[1]);converted=binary(sources[2]);loaded=binary(sources[3]);final=binary(sources[4])
    np.testing.assert_array_equal(converted['GlobalIndex'],np.arange(5233))
    np.testing.assert_array_equal(np.column_stack([original['x'],original['y']]),mesh_coordinates(MESH))
    for name in ['x','y']:
        np.testing.assert_array_equal(original[name],loaded[name]);np.testing.assert_array_equal(final[name],loaded[name])
    rho=original['Density'];mom=np.column_stack([original['Momentum_x'],original['Momentum_y']])
    pressure=(1.4-1)*(original['Energy']-.5*np.sum(mom**2,axis=1)/rho)
    np.testing.assert_array_equal(loaded['Density'],rho);np.testing.assert_allclose(loaded['Pressure'],pressure,rtol=4e-16,atol=0)
    np.testing.assert_array_equal(converted['Pressure'],loaded['Pressure'])
    assert np.min(final['Density'])>0 and np.min(final['Pressure'])>0
    geometry=geometric_faces();native=align(geometry,surface(sources[6]));primitive=align(geometry,surface(sources[7]))
    # Same surface quadrature, matching geometry and owner: no interpolation.
    for x,y in zip(native,primitive):
        for k in ['marker','node','x','y','nx','ny']:assert x[k]==y[k]
    rows=[]
    for f,n,w in zip(geometry,native,primitive):
        # Export the actual common sample coordinates, after the independent
        # geometric check (alternative arithmetic differs by <= roundoff).
        f=dict(f,x=n['x'],y=n['y'])
        rows.append(dict(**f,cp_native=n['Cp'],cp_primitive=w['Cp'],delta_cp=(w['p']-n['p'])/QINF,
            p_native=n['p'],p_primitive=w['p'],fx_native=(n['p']-PINF)*n['nx'],fx_primitive=(w['p']-PINF)*w['nx'],
            delta_fx=(w['p']-n['p'])*n['nx'],delta_cd=(w['p']-n['p'])*n['nx']/(QINF*AREA)))
    out.mkdir(parents=True,exist_ok=False)
    write(out/'cp_comparison.csv',rows)
    metrics=[]
    for side in ['all','upper','lower']:
        for group in REGIONS:
            selected=[r for r in rows if (side=='all' or r['side']==side) and (group=='all' or r['region']==group)]
            d=np.array([r['delta_cp'] for r in selected]);length=np.array([r['length'] for r in selected]);i=int(np.argmax(abs(d)))
            metrics.append(dict(side=side,region=group,samples=len(d),length=sum(length),mean=float(np.mean(d)),
                rms=float(np.sqrt(np.mean(d*d))),maximum=float(np.max(abs(d))),x_at_max=selected[i]['x'],
                y_at_max=selected[i]['y'],length_weighted_mean=float(np.dot(length,d)/sum(length)),
                length_weighted_rms=float(np.sqrt(np.dot(length,d*d)/sum(length)))))
    write(out/'cp_metrics.csv',metrics)
    total_delta=math.fsum(r['delta_cd'] for r in rows)
    decomposition=[]
    for group in REGIONS:
        selected=[r for r in rows if group=='all' or r['region']==group]
        fn=math.fsum(r['fx_native'] for r in selected);fw=math.fsum(r['fx_primitive'] for r in selected)
        df=math.fsum(r['delta_fx'] for r in selected)
        decomposition.append(dict(region=group,samples=len(selected),fx_native=fn,fx_primitive=fw,delta_fx=df,
            delta_cd=df/(QINF*AREA),signed_percent=100*df/(QINF*AREA*total_delta)))
    write(out/'drag_decomposition.csv',decomposition)
    np.testing.assert_allclose(math.fsum(r['delta_cd'] for r in decomposition[1:]),total_delta,rtol=1e-14,atol=1e-16)
    archived=history(folder/'history.csv');zero=archived[0];last=archived[-1];assert last['Time_Iter']==9000
    fxn,fxw=[math.fsum(r[k] for r in rows) for k in ['fx_native','fx_primitive']]
    fyn=math.fsum((r['p_native']-PINF)*r['ny'] for r in rows);fyw=math.fsum((r['p_primitive']-PINF)*r['ny'] for r in rows)
    for computed,expected in [(fxn/QINF,zero['CD']),(fxw/QINF,last['CD']),(fyn/QINF,zero['CL']),(fyw/QINF,last['CL'])]:
        np.testing.assert_allclose(computed,expected,rtol=0,atol=3e-14)
    # Independent native convention: nodal pressure times the sum of its
    # geometric half-face normals equals native vertex-normal integration.
    nodal_fx=math.fsum((pressure[r['node']]-PINF)*r['nx'] for r in rows)
    nodal_fy=math.fsum((pressure[r['node']]-PINF)*r['ny'] for r in rows)
    native_last=history(sources[8])[-1]
    np.testing.assert_allclose([nodal_fx/QINF,nodal_fy/QINF],[native_last['CD'],native_last['CL']],rtol=0,atol=3e-14)
    write(out/'forces.csv',[dict(field='native_common_reconstructed',Fx=fxn,Fy=fyn,CD=fxn/QINF,CL=fyn/QINF),
                            dict(field='primitive_9000_reconstructed',Fx=fxw,Fy=fyw,CD=fxw/QINF,CL=fyw/QINF),
                            dict(field='native_nodal_independent_check',Fx=nodal_fx,Fy=nodal_fy,CD=nodal_fx/QINF,CL=nodal_fy/QINF)])
    tip_segments={r['segment'] for r in rows if r['node']==199}
    tip=[dict(**r,owner_is_te=r['node']==199) for r in rows if r['segment'] in tip_segments]
    assert len(tip)==4 and sum(r['owner_is_te'] for r in tip)==2
    write(out/'te_adjacent_half_faces.csv',tip)
    # Fine TE bins keep cancellations visible, without modifying the coarse partition.
    fine=[]
    for left,right in [(.9,.98),(.98,.99),(.99,.999),(.999,1.)]:
        selected=[r for r in rows if left<r['x']<=right]
        fine.append(dict(x_min=left,x_max=right,samples=len(selected),delta_cd=math.fsum(r['delta_cd'] for r in selected)))
    write(out/'trailing_edge_bins.csv',fine)
    assert all(hashlib.sha256(path.read_bytes()).hexdigest()==hashes[str(path.relative_to(ROOT))] for path in sources)
    (out/'provenance.json').write_text(json.dumps(dict(kind='diagnostic_only_not_stationary_certification',
        sources_sha256=hashes,pressure_precision='17-significant-digit CSV from Float64 solver fields',
        interpolation='none',p_inf=PINF,q_inf=QINF,ref_area=AREA,samples=400,
        upper=sum(r['side']=='upper' for r in rows),lower=sum(r['side']=='lower' for r in rows),
        native_cd_common=fxn/QINF,primitive_cd=fxw/QINF,delta_cd=total_delta,
        cd_reproduction_error=max(abs(fxn/QINF-zero['CD']),abs(fxw/QINF-last['CD'])),
        spatial_rhs_final=last['RHS_CURRENT'],te_rhs_final=last['RHS_TE']),indent=2))
    plots(out,rows,decomposition)
    print(json.dumps(dict(delta_cd=total_delta,regions=decomposition,metrics=metrics),indent=2))


def plots(out,rows,decomposition):
    os.environ.setdefault('MPLCONFIGDIR','/private/tmp/primitive-phase15-matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'savefig.dpi':300,'axes.grid':True,'grid.alpha':.25,'pdf.fonttype':42})
    # Color encodes field, line/marker encodes surface; each branch is sorted alone.
    def curves(ax,limits,markers=False):
        for field,color in [('native','#1764ab'),('primitive','#cc3b26')]:
            for side,style,marker in [('upper','-','o'),('lower','--','s')]:
                selected=sorted([r for r in rows if r['side']==side and limits[0]<=r['x']<=limits[1]],key=lambda r:r['x'])
                ax.plot([r['x'] for r in selected],[r['cp_'+field] for r in selected],color=color,linestyle=style,
                    marker=marker if markers else None,markersize=3,linewidth=1.2,
                    label=('Native field, common quadrature' if field=='native' else 'Primitive, 15D / 9000 (not certified steady)')+' — '+side)
        ax.set(xlabel='$x/c$',ylabel='$C_p$',xlim=limits);ax.invert_yaxis()
    fig,ax=plt.subplots(figsize=(10,6));curves(ax,(0,1));ax.legend(fontsize=8)
    ax.set_title('NACA 0012 — M=0.2, AOA=0°\nCommon reconstruction/quadrature; primitive not certified stationary')
    fig.tight_layout();fig.savefig(out/'cp_native_vs_primitive.png');fig.savefig(out/'cp_native_vs_primitive.pdf');plt.close(fig)
    for name,limits in [('cp_trailing_edge_zoom',(.9,1)),('cp_trailing_edge_detail',(.98,1)),('cp_leading_edge_zoom',(0,.1))]:
        fig,ax=plt.subplots(figsize=(9,5));curves(ax,limits,True);ax.legend(fontsize=7)
        ax.set_title('Native vs primitive 15D / 9000 — diagnostic, not certified steady')
        fig.tight_layout();fig.savefig(out/(name+'.png'));fig.savefig(out/(name+'.pdf'));plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    for ax,limits in zip(axes,[(0,1),(.9,1)]):
        for side,color,style in [('upper','#1764ab','-'),('lower','#cc3b26','--')]:
            selected=sorted([r for r in rows if r['side']==side and limits[0]<=r['x']<=limits[1]],key=lambda r:r['x'])
            ax.plot([r['x'] for r in selected],[r['delta_cp'] for r in selected],style,color=color,marker='.',markersize=3,label=side)
        ax.axhline(0,color='black',linewidth=.6);ax.set(xlabel='$x/c$',ylabel='$\\Delta C_p$ (primitive − native)',xlim=limits);ax.legend()
    fig.suptitle('Common points; no interpolation — primitive not certified stationary');fig.tight_layout()
    fig.savefig(out/'delta_cp.png');fig.savefig(out/'delta_cp.pdf');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5));values=[r['delta_cd'] for r in decomposition[1:]]
    bars=ax.bar(['LE: x<0.10','Central: 0.10≤x≤0.90','TE: x>0.90'],values,color=['#1764ab','#777777','#cc3b26'])
    for bar,val in zip(bars,values):ax.annotate(f'{val:+.6g}',(bar.get_x()+bar.get_width()/2,val),ha='center',va='bottom' if val>=0 else 'top')
    ax.axhline(0,color='black',linewidth=.8);ax.set(ylabel='Contribution to $\\Delta C_D$',title='Axial-force difference, common half-face quadrature')
    fig.tight_layout();fig.savefig(out/'drag_difference_by_region.png');fig.savefig(out/'drag_difference_by_region.pdf');plt.close(fig)
    ordered=sorted(rows,key=lambda r:r['x']);cumulative=np.cumsum([r['delta_cd'] for r in ordered])
    fig,ax=plt.subplots(figsize=(8,5));ax.step([r['x'] for r in ordered],cumulative,where='post')
    ax.set(xlabel='$x/c$',ylabel='Cumulative $\\Delta C_D$',title='Sum of discrete half-face contributions, both surfaces')
    fig.tight_layout();fig.savefig(out/'drag_difference_cumulative.png');plt.close(fig)
    write(out/'drag_cumulative.csv',[dict(x=r['x'],side=r['side'],node=r['node'],cumulative_delta_cd=float(c)) for r,c in zip(ordered,cumulative)])


if __name__=='__main__':main()
