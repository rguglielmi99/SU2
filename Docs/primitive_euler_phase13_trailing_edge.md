# FASE 13 — Topologia e vincoli al trailing edge

## Audit e intervento

Modificati `CPrimitiveEulerSolver.cpp`, solo `BuildSlipConstraints`, e
`primitive_euler_boundary_tests.cpp`, quattro test `[primitive-trailing-edge]`.
Nuovo file: questo rapporto. Nessuna modifica a geometria nativa, normali,
kernel, riflessione, fluttuazioni, gradienti, MUSCL, farfield o FE.

`CAggregatedGeometry::BoundaryFace` conserva nodo, marker, punto e normale
uscente integrata della mezza faccia; la connettività è basata sugli **ID**,
non sull'uguaglianza delle coordinate. La normale nodale FASE 10 viene
costruita per coppia (nodo,marker). Gli ID duplicati non sono fusi.

Il modello precedente assumeva implicitamente liscio ogni marker Euler.
Due normali quasi opposte potevano produrre una risultante non nulla ma
parallela alla corda: la proiezione eliminava tutta la velocità longitudinale.
Separando i marker si otteneva invece rango due sullo stesso valore nodale,
ancora con velocità nulla. Sono meccanismi distinti dal caso di risultante nulla.

## Politica minima implementata

1. **Stesso marker dichiarato liscio:** la media FASE 10 resta invariata
   per normali coorientate, n1·n2>0 entro roundoff. Sono rifiutate le pieghe
   non coorientate, con messaggio `fold/trailing edge` e ID del nodo.
   È una restrizione conservativa di ammissibilità della media, **non** un
   classificatore universale di spigoli o una soglia fisica di Kutta.
   Per due facce coorientate ||S1+S2||/(|S1|+|S2|)>=1/sqrt(2): la risultante
   è ben condizionata. Curve troppo grossolanamente campionate possono essere
   rifiutate e richiedere raffinamento; non viene introdotto un angolo regolabile.
2. **Marker distinti, vincoli indipendenti:** si sommano gli angoli dei
   triangoli fluidi incidenti tramite `nodes->GetElems`. Se il settore fluido
   supera pi, si rifiuta il vertice reentrante con diagnostica esplicita:
   rango due imporrebbe stagnazione. Rimangono supportati gli angoli convessi
   del dominio fluido già verificati nella scatola, con rango due intenzionale.
3. **Normali opposte ma collineari su marker distinti:** rango uno; il
   vincolo già esistente consente la tangente comune. Nessuna proiezione speciale.
4. Risultante nulla, non rappresentabile o cancellata entro roundoff:
   rimane il rifiuto `degenerate resultant wall normal`. Rimangono anche i
   precedenti controlli non-manifold e sulle intersezioni di marker non rettilinei.

Le tolleranze 256 epsilon identificano il roundoff, non una classe fisica
di trailing edge. Gli spigoli autentici devono essere dichiarati tramite
marker distinti: uno spigolo poco marcato nello stesso marker resta
indistinguibile da un campionamento di curva senza ulteriori metadati.
Non si pretende di rilevare automaticamente ogni spigolo geometrico.

Il rifiuto avviene sia prima della proiezione FREESTREAM sia prima del
residuo. I test verificano che nessuno stato venga modificato prima dell'errore.

## Misure geometriche

S sono normali di mezze facce; r=||S1+S2||/(|S1|+|S2|).
R è il rango delle **due normali geometriche**, distinto dal rango uno
del modello mediato per parete liscia. Le velocità candidate riportate
ricostruiscono algebricamente il comportamento precedente; sui casi rifiutati
non viene più applicata alcuna proiezione produttiva.

| Caso | S1, S2 | Somma | Angolo normali | r | R / esito |
|---|---|---|---:|---:|---|
| Curva FASE 10 | (.1,-.5),(.05,-.5) | (.15,-1) | 5.60° | circa .999 | 2 / media rango 1 invariata |
| TE acuto | (-.06,.3),(-.06,-.3) | (-.12,0) | 157.38° | .196116 | 2 / rifiuto |
| Quasi cuspide | (-5e-9,.3),(-5e-9,-.3) | (-1e-8,0) | 179.999998° | 1.667e-8 | 2 / rifiuto |
| Cuspide ideale, marker distinti | (0,.3),(0,-.3) | (0,0) | 180° | 0 | 1 / tangente comune |
| Chiusura spessa, angolo inferiore | (-.03,.5),(-.06,0) | (-.09,.5) | 86.57° | .905752 | 2 / settore reentrante rifiutato |

Con freestream (.4,0), nel TE acuto la media candidata dà (0,0), pur con
risultante ben rappresentabile. Prima della proiezione le componenti locali
sono u_n=-.0784465 e u_t=±.392232; dopo sarebbero entrambe nulle.
La quasi cuspide ha risultante molto piccola rispetto alle facce ma ancora
superiore alla tolleranza di degenerazione: viene rifiutata per topologia
non supportata, non perché la normalizzazione sia numericamente impossibile.

La chiusura spessa ha due vertici separati di .12 e un segmento di base:
dividerla in tre marker espone due vincoli indipendenti per ciascun angolo,
non elimina il problema. Il test non usa una media per occultarlo.

Due nodi coincidenti con ID distinti, in due triangoli fluidi separati,
conservano invece vincoli indipendenti per ID: la proiezione dà
(.396040,-.039604) e (.396040,+.039604), con velocità tangenziale non nulla.
È un test topologico, non una mesh di profilo pronta all'uso.

## Cuspide e verifica dinamica

Una cuspide ideale può avere tangente limite comune; un vertice ad angolo
finito con due normali indipendenti e un solo vettore velocità soddisfa
entrambi gli slip soltanto con u=v=0. Questa algebra non autorizza una
stagnazione universale al trailing edge.

Un dominio anulare triangolato con fessura di spessore nullo, lati con ID
intermedi distinti e tip condiviso, rappresenta la cuspide ideale nel test:
le due normali sono opposte, rango uno. Tre passi FE preservano lo stato
uniforme tangenziale (.4,0), positivo, con RHS nullo entro roundoff.
Anche il caso di ID coincidenti distinti esegue tre passi: stati/RHS finiti,
rho/p positive, slip preservato e velocità tangenziale >.3.
Nessuna condizione di Kutta o scelta della circolazione è implementata.

## Mesh NACA locale e raccomandazione

Esaminata senza simulazione `QuickStart/mesh_NACA0012_inv.su2`, stesso nome
referenziato dal cfg NACA nei TestCases: 10216 triangoli, marker airfoil
200 segmenti, farfield 50. Il TE è il solo ID **199**, (1,0), vicino agli ID
0/198 con coordinate (.99975001812,±3.632896519e-5). Nessun ID duplicato sul
profilo. Normali unitarie incidenti (-.143816,±.989604), angolo 163.4626°,
risultante lungo -x. Il solo nodo airfoil che viola la coorientazione è 199.
Questa mesh è quindi **intenzionalmente non supportata** dalla nuova politica.
SHA256 del file esaminato:
`9094b51c2628d3bb4c865d774b2308dbb59a213ad19c007cfcd7d57e7aaeebeb`.

Per il NACA originale raccomando di conservare questa topologia come caso
di sviluppo, senza duplicare ID o separare marker come espediente. Servirà
una formulazione TE esplicita; eventuali DOF distinti con taglio di scia
richiedono connettività e accoppiamento ulteriori, non solo coordinate uguali.
Per uno studio preliminare già nel perimetro geometrico, usare una chiusura
arrotondata ben risolta e dichiarare che modifica il profilo originale;
una base piatta con angoli vivi non è automaticamente una soluzione.

## Riproducibilità ed esiti

```sh
./ninja -C build SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_tests.py --output-dir /private/tmp/primitive-phase13-build
/private/tmp/primitive-phase13-build/test_driver '[primitive-trailing-edge]'
/private/tmp/primitive-phase13-build/test_driver '[primitive-boundary],[primitive-wall-gradients],[primitive-geometry],Rotating disk measures*,Global balance chain*,Constant primitive fields*,Affine NODE*,Forward Euler applies*,Rejected time steps*'
git diff --check
```

Build riuscita; **4 nuovi casi/574 asserzioni**, **27 regressioni/22106
asserzioni**, tutti passati. Stampe TE_AUDIT/TE_COMPONENT riproducono le
misure. Nessuna campagna completa MMS, vortice o conservazione globale.

**A:** curve lisce ammissibili, angoli convessi fluidi dichiarati, cuspide
ideale testata e ID distinti nel caso topologico ristretto.
**B:** pieghe non coorientate, TE acuti condivisi, angoli reentranti con
vincoli indipendenti, degenerazioni.
**C:** TE aerodinamico NACA, wake coupling e circolazione richiedono sviluppo.
Nessuna validazione aerodinamica o modifica Git dello stato.
