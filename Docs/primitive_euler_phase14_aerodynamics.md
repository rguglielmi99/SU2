# FASE 14 — Coefficienti aerodinamici e output

## Perimetro e audit nativo

Branch `feature/primitive-path-solver`, HEAD iniziale `04edb84844`, working tree
inizialmente pulito. Nessuna operazione Git di scrittura. Nessuna campagna NACA
stazionaria, MMS, vortice o bilanci globali.

- [CFVMFlowSolverBase::Pressure_Forces](../SU2_CFD/include/solvers/CFVMFlowSolverBase.inl#L1852)
  legge pressione **nodale**, `GetMarker_All_Monitoring`, normali native ai
  boundary vertex e accumula `-(p-p_inf)*Normal/AeroCoeffForceRef`.
  La geometria primitiva usa invece `S=-Normal` dopo somma delle mezze facce
  incidenti allo stesso nodo/marker: stesso segno fisico della forza sul corpo.
- [CEulerSolver::SetReferenceValues](../SU2_CFD/src/solvers/CEulerSolver.cpp#L1434)
  definisce `DynamicPressureRef=.5*rho_inf*|V_inf|²` e
  `AeroCoeffForceRef=DynamicPressureRef*REF_AREA`. Il fallback nativo a 1
  per pressione dinamica piccola **non è riprodotto**: qui q non positivo
  viene rifiutato.
- [CFVMFlowSolverBase](../SU2_CFD/include/solvers/CFVMFlowSolverBase.hpp#L137)
  conserva `InvCoeff`, `SurfaceInvCoeff`, `AllBoundInvCoeff` e
  `CPressure[marker][vertex]`. Queste strutture non sono usate dal primitivo.
- [CFlowOutput](../SU2_CFD/src/output/CFlowOutput.cpp#L1814) registra
  `LIFT`/`DRAG`, etichette CL/CD, gruppo `AERO_COEFF`, e legge i getter del
  solver nativo; [CFlowCompOutput](../SU2_CFD/src/output/CFlowCompOutput.cpp#L346)
  esporta Cp puntuale nodale tramite il normale percorso dei campi.
- [CConfig](../Common/src/CConfig.cpp#L1540): `REF_AREA` normalizza le forze;
  `REF_LENGTH` normalizza i **momenti**, non va moltiplicato nuovamente
  nel denominatore di CL/CD. Entrambi hanno default 1.

Riutilizzati registrazione history, selezione nominale `MARKER_MONITORING`,
AOA in gradi, REF_AREA e calendario `COutput::SetResultFiles` /
`WriteAdditionalFiles`. Nessun percorso EULER/NAVIER_STOKES modificato.

## Calcolo implementato

`CPrimitiveEulerSolver::ComputeAerodynamics` restituisce forza, CL/CD e
campioni superficiali; richiama **PrepareReconstruction esistente** sullo stato
corrente, poi `ReconstructAt(node,q)` per ogni mezza faccia selezionata.
Questo evita gradienti dello stadio precedente dopo un aggiornamento FE.
Non modifica W, tempo o LinSysRes; aggiorna soltanto le cache di ricostruzione
e i relativi limiter/fattori di positività. La diagnostica MIN_THETA nell'output
monitorato riguarda quindi questa ricostruzione corrente, mentre RHS rimane
quello dello stadio accettato.

Con S uscente dal fluido:

    F = sum_b (p_b-p_inf) S_b
    q_inf = .5 rho_inf (u_inf²+v_inf²)
    Cp_b = (p_b-p_inf)/q_inf
    CD = (Fx cos(alpha)+Fy sin(alpha))/(q_inf REF_AREA)
    CL = (-Fx sin(alpha)+Fy cos(alpha))/(q_inf REF_AREA)

In 2D REF_AREA è corda × apertura unitaria: per NACA di corda 1, impostare
`REF_AREA=1`, denominatore q_inf. Dati fisici SI e scale DIMENSIONAL:
freestream ottenuto dallo **stesso adapter delle BC primitive** (rho da
DENSITY_FS oppure p/(RT), velocità da Mach, gamma e AOA). Non viene dedotto
da `PRIMITIVE_INITIAL_STATE`, che può rappresentare un altro stato.
Rifiutati q<=0, riferimenti non finiti, REF_AREA<=0 e overflow dei coefficienti;
nessun calcolo automatico di REF_AREA=0. CM e coefficienti viscosi esclusi.

La selezione deve contenere solo MARKER_EULER esistenti: CUSTOM/FAR e duplicati
sono rifiutati, non integrati silenziosamente. Più marker condividono nodi ma
non campioni: ogni mezza faccia geometrica contribuisce una volta.
Senza MARKER_MONITORING l'API restituisce `enabled=false`, superficie vuota;
output aero disabilitato e nessuna normalizzazione richiesta.

La routine non applica vincoli di slip: consente verifiche geometriche di
pressione anche sulla NACA originale, senza autorizzare il TE debole nel
driver. Non certifica che lo stato assegnato sia una soluzione delle BC.

## Uso e output ordinario

Configurazione esistente, nessuna nuova opzione:

```ini
MARKER_MONITORING=(airfoil)
REF_AREA=1
HISTORY_OUTPUT=(ITER,TIME_DOMAIN,PRIMITIVE_RHS,PRIMITIVE_STATE,AERO_COEFF)
```

Se HISTORY_OUTPUT non è specificato, AERO_COEFF è incluso automaticamente
quando vi sono marker monitorati. Si possono richiedere anche LIFT/DRAG
individualmente. CSV superficiale:
`<SURFACE_FILENAME>_primitive_pressure_<Time_Iter>.csv`, alle scadenze dei file
di risultato, una volta per indice temporale anche con più formati richiesti.
Nessuna scrittura durante valutazioni del residuo. Campioni di mezza faccia,
non Cp nodali interpolati tramite il writer di superficie nativo.

Esempio reale del piccolo caso CLI, indice 5:

```csv
marker,node,x,y,p,Cp,nx,ny
"wall",0,0.25,0.05,0.9998605905819381,-0.0013940941806189764,0.1,-0.5
```

History a T=.0005: CL=.0033974476006990001,
CD=.001113475496358475. Sono coefficienti del piccolo **bordo aperto**
poligonale del test CLI, non risultati NACA. La somma indipendente del CSV
riproduce entrambi entro 1e-12.

## Verifiche analitiche e regressioni

Corpo quadrato unitario circondato da fluido, due marker di corpo e uno esterno
escluso; q=.08, REF_AREA=2, alpha=30°. WLS affine, limiter NONE.

| Campo p | Fx | Fy | CL | CD |
|---|---:|---:|---:|---:|
| 2 uniforme | 0 | 0 | 0 | 0 |
| 2+.2x-.1y | -.2 | .1 | 1.166265877365274 | -.770031754730548 |

Forze verificate entro 2e-13, coefficienti entro 3e-12; uniforme Cp=12.5.
Controllati pressione al punto di quadratura, somma manuale, copertura univoca,
normali verso il corpo, chiusura, inversione dell'ordine dei segmenti e
connettività, aggiornamento da uniforme ad affine senza preparazione esterna,
LinSysRes/tempo invariati. Testati assenza selezione, marker non-wall,
duplicati, REF_AREA zero/negativa e Mach zero.

NACA originale: 200 segmenti/400 campioni airfoil, due campioni distinti al
tip globale 199, pressione uniforme 2, Cp=12.5, |Fx| e |Fy|<2e-13.
Controllo aggiunto al test TE preesistente **prima** dell'autorizzazione debole;
poi ripristinato lo stato originale e rieseguiti invariati i cinque passi.

```sh
./ninja -C build SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/primitive-phase14-build
/private/tmp/primitive-phase14-build/test_driver '[primitive-aerodynamics]'
/private/tmp/primitive-phase14-build/test_driver '[primitive-geometry],[primitive-boundary],[primitive-trailing-edge],[primitive-wall-gradients],Constant primitive fields*,Affine NODE*,Forward Euler applies*,Rejected time steps*'
/private/tmp/primitive-phase14-build/test_driver '[primitive-weak-te]'
python3 UnitTests/check_primitive_euler_polygon_cli.py --executable build/SU2_CFD/src/SU2_CFD --output-dir /private/tmp/primitive-phase14-cli
git diff --check
```

Build riuscita; rispettivamente **1 caso/218 asserzioni**, **29/18002**,
**1/315673**, tutti passati. CLI: cinque passi, quattro snapshot e quattro
tabelle; positività/slip conservati, history confrontata con i CSV.

File modificati: header/source di CPrimitiveEulerSolver e CPrimitiveEulerOutput;
primitive_euler_boundary_tests.cpp; check_primitive_euler_polygon_cli.py;
questo rapporto. Nessuna modifica a kernel, geometria, BC, ricostruzione o FE.
Il TE debole rimane sperimentale solo API. Calcolo dei coefficienti verificato;
accuratezza fluidodinamica, portanza/resistenza convergenti e stazionarietà
su profilo restano da validare nella FASE 15. Nessun commit o push.
