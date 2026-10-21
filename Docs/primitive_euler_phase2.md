# Euler 2D primitivo — FASE 2: AGGREGATED e NODE

Questa fase prepara le metriche geometriche per il futuro solver. Non aggiunge
un solver eseguibile, opzioni di configurazione o flussi al bordo e non modifica
la geometria o il solver Euler conservativo SU2.

## API

Header: `SU2_CFD/include/geometry/primitive_euler_geometry.hpp`.

```cpp
// Dopo il preprocessing geometrico standard, inclusi SetEdges e SetControlVolume.
PrimitiveEuler2D::CAggregatedGeometry metrics(*geometry, *config);

const auto& face = metrics.GetFace(iEdge);
// face.first, face.second: indici LOCALI degli estremi, come in SU2.
// face.point: Qeff=(xi+xj)/2.
// face.normal: normale integrata Sij di SU2, da first verso second.
const auto d = PrimitiveEuler2D::ComputeFluctuations(Wi, Wj, face.normal, gamma);

const auto& cell = metrics.GetCell(iPoint);
// cell.volume: volume duale originale SU2.
// cell.offset: centroide geometrico del duale meno coordinata del nodo.
```

La normale non è unitaria: `ComputeFluctuations` applica già la misura della
faccia. Non moltiplicare nuovamente per la sua norma. Qeff è il punto medio
dell'edge primale, non il centroide della faccia duale spezzata.

Gli stati restano nodali, `Wi = W(xi)`. I centroidi servono a rendere possibile
la chiusura MEAN_COEFF della specifica, anche in NODE:

```
Wmean = Wi + Glim_i * cell.offset
```

Questa formula sarà usata nel residuale di volume di una fase successiva.
Il calcolo dei centroidi non implementa DUAL_MEAN, e non costituisce una prova
dell'ordine due dello schema.

## Costruzione e controlli

La classe è una fotografia immutabile, costruita senza scrivere in `CGeometry`.
Copre solo mesh statiche, seriali, triangolari 2D, sul livello più fine e senza
periodicità. Rifiuta mesh con halo, livelli multigrid, altre dimensioni, elementi
non triangolari, geometria dinamica e configurazioni adjoint. Non è validata per
tipi AD. Ricostruire la fotografia se cambiano coordinate o topologia.

Per ogni triangolo si usano gli stessi sottotriangoli di SU2, con vertici
`(nodo, punto medio edge, centro triangolo)`. La loro area e il primo momento
permettono di ricavare il centroide duale, inclusi i nodi di bordo. I momenti
sono accumulati relativamente al nodo per ridurre cancellazioni dovute
all'origine delle coordinate. Il centro del triangolo primale non viene
confuso con il centroide duale.

Le aree ricostruite vengono confrontate con i volumi SU2; la somma dei
contributi geometrici delle normali viene confrontata con Sij. La tolleranza è
relativa, `1e-11`, rispettivamente all'area duale e alla somma delle lunghezze
dei contributi per edge: non introduce una soglia assoluta dipendente dalle
unità della mesh. Si controllano inoltre orientazione, incidenza degli edge,
volumi positivi e valori finiti. Un'incoerenza genera `std::invalid_argument`.
I volumi e le normali effettivamente esposti restano quelli originali SU2.

Non viene memorizzata né esposta alcuna lista di segmenti, punti Qm o modalità
di quadratura SEGMENTED. La scomposizione locale serve esclusivamente per i
momenti geometrici e il controllo delle metriche aggregate.

## Verifica

`UnitTests/SU2_CFD/primitive_euler_geometry_tests.cpp`, tag `[primitive-geometry]`,
crea mesh SU2 temporanee e usa `CPhysicalGeometry`, il reader e il preprocessing
standard, compresi `SetControlVolume` e `SetBoundControlVolume`.

Sono verificati:

- triangolo singolo: volumi, normali integrate e punto Qeff;
- ventagli triangolari regolari e distorti: aree e centroidi contro la formula
  indipendente di shoelace sui poligoni duali;
- area totale e primi momenti dell'intero dominio;
- media esatta di un campo affine sui duali, inclusi quelli di bordo;
- chiusura dei duali interni e di bordo, includendo i vertici d'angolo;
- rinumerazione dei nodi, inversione delle normali e corrispondenti fluttuazioni;
- rifiuto di metriche alterate, triangoli invertiti e modalità non supportate.

Le normali memorizzate nei vertici di bordo SU2 risultano entranti; il controllo
di chiusura le sottrae alla somma delle normali interne uscenti. Questo controllo
geometrico non implementa né valida un trattamento path-conservative del bordo.

Il file è registrato in `UnitTests/meson.build`. Con il driver SU2 abilitato:

```sh
./meson.py configure build -Denable-tests=true
./ninja -C build UnitTests/test_driver
build/UnitTests/test_driver '[primitive-geometry],[primitive-euler]'
```

Per questa consegna il driver dedicato in `/private/tmp/su2-phase2-tests` è
stato compilato con il compilatore, le macro e i flag della build esistente
(`build/compile_commands.json`), e collegato alle sue librerie SU2/MPI. Non è
stata riconfigurata la build principale, che conserva `enable-tests=false`.

```sh
python3 /private/tmp/su2-build-phase2-tests.py
/private/tmp/su2-phase2-tests '[primitive-geometry],[primitive-euler]'
```

Gli artefatti temporanei non sono parte del repository. Non sono stati eseguiti
l'intera suite SU2, studi MMS o simulazioni del futuro solver. SEGMENTED,
DUAL_MEAN, mesh miste, MPI distribuito e griglie mobili rimangono rinviati.

Esito della consegna: **20 test superati, 137894 asserzioni**, di cui 5 test
geometrici nuovi e 15 test del kernel. Compilazione Apple Clang 15, C++17,
con i flag della build SU2 esistente. Sono comparsi avvisi non bloccanti già
legati all'ambiente/infrastruttura: pragma GCC non riconosciuto da Clang,
librerie duplicate nella link line SU2 e bind TCP MPI non consentito nel
sandbox. L'esecuzione seriale e i test sono terminati con codice zero.

## Verifiche geometriche aggiuntive prima della FASE 3

La suite aggiunge una triangolazione deterministica del rettangolo `[0,3] x
[0,2]`: 9 nodi, 18 edge, 10 triangoli CCW di aree e forme non uniformi. Il
contorno ha quattro angoli e due nodi su lati rettilinei, `(1.1,0)` e `(1.7,2)`.
I tre nodi strettamente interni sono `(0.8,0.8)`, `(2.2,1.1)` e `(1.4,0.55)`,
con valenze rispettivamente 6, 6 e 3. Coordinate e connettività sono esplicite
in `irregularXY` e `irregularTriangles` nel file dei test. Vengono controllate
l'orientazione iniziale positiva, le valenze, la variazione delle aree e le
dimensioni della connettività. Non è una semplice deformazione del ventaglio
precedente.

**Riferimenti indipendenti.** Per ciascun triangolo, i contributi ai duali
sono quadrilateri `(nodo, midpoint successivo, centro triangolo, midpoint
precedente)`, integrati con shoelace. Area e centroide sono calcolati in
coordinate relative al primo vertice, con accumulo `long double`. Per ogni
duale si accumulano momenti relativi al suo nodo, evitando prodotti fra area
e grandi coordinate assolute. Si confrontano direttamente gli offset,
senza ricostruire due centroidi assoluti per poi sottrarli. Per il dominio,
un secondo riferimento usa esclusivamente il poligono esterno: area 6 e
primi momenti assoluti 9 e 6. Si verificano anche i momenti rispetto a
un'origine locale, per rendere significativo il controllo dopo traslazione.

**Controlli.** Tutti i nodi sono verificati per volume, offset e media
geometrica di due campi affini, con gradienti `(0.3,-0.2)` e `(-0.7,0.6)`.
Si controlla anche la correzione affine isolata dal valore nodale, per non
mascherare errori piccoli con un grande valore di fondo. La chiusura viene
verificata separatamente per i tre interni, i due nodi di bordo rettilineo e
i quattro angoli. Le normali native di bordo sono confrontate direttamente
con le seminomali entranti costruite dal contorno, oltre al controllo di
chiusura.

**Convenzioni confermate nel sorgente SU2.**

- `Common/src/geometry/CGeometry.cpp`, `SetEdges` (riga 1650): chiama
  `edges->SetNodes(iEdge,iPoint,jPoint)` solo dentro `if (iPoint < jPoint)`.
  L'ordine locale crescente è quindi garantito per la geometria supportata
  dopo il preprocessing standard. Non è una convenzione sugli ID globali.
- `Common/src/geometry/CPhysicalGeometry.cpp`, `SetControlVolume` (riga
  6683): `change_face_orientation = (face_iPoint > face_jPoint)` orienta i
  contributi coerentemente con quell'ordine. `CEdge::SetNodes_Coord` 2D
  (`Common/src/geometry/dual_grid/CEdge.cpp`, riga 84) accumula `(dy,-dx)`.
- `SetBoundControlVolume` (riga 6836) e `CVertex::SetNodes_Coord` 2D
  (`Common/src/geometry/dual_grid/CVertex.cpp`, riga 48) producono le
  normali entranti ai vertici dopo i controlli di orientazione del bordo.
  Per la chiusura si somma `+Sij` al primo nodo e `-Sij` al secondo e si
  **sottrae** la normale nativa di bordo. Il test confronta questa formula
  con seminomali analitiche `(-dy,dx)/2` sui lati CCW del dominio.

**Scale e tolleranze.** Il controllo interno `RequireClose` rimane invariato:
`1e-11` relativo all'area duale o alla somma delle lunghezze geometriche.
Il helper `Check` verifica anche finitezza del riferimento e positività della
scala; nei nuovi confronti le scale dimensionali sono lunghezza, lunghezza
al quadrato e lunghezza al cubo per offset/normali, aree e momenti. La soglia
di test resta `2e-12 * max(scala, |riferimento|)`, senza allargamenti per far
passare i test. Sono verificati il riscalamento uniforme `1e-3` e il
riscalamento 8 con traslazione `(128,-256)`: volumi proporzionali al quadrato,
normali e offset proporzionali alla scala, Qeff trasformato coerentemente.
Un test distinto verifica il solo riferimento shoelace con origine
`(1e9,-1e9)` su un rettangolo rappresentabile esattamente; non pretende che
l'intera geometria SU2 sia robusta a quel rapporto fra origine e dimensioni.
Restano possibili perdite di precisione nei centri assoluti calcolati da
SU2 quando tale rapporto è molto grande, e nelle mesh quasi degeneri.

**Controlli negativi sicuri.** Dopo il preprocessing di una mesh valida, si
alterano le coordinate per costruire un triangolo collineare o una coordinata
NaN; separatamente si sostituiscono elementi con connettività già esistente
per avere tre triangoli incidenti sullo stesso edge. Si controllano i messaggi
specifici dei relativi guard. È aggiunta anche un'alterazione del modulo della
normale nativa del 25%, oltre all'inversione già testata. Il reader SU2 non
riceve dati degeneri o non finiti e non vengono alterate allocazioni o indici
fuori intervallo.

**Esito effettivo.** Tutti i 20 test precedenti restano presenti. La suite
pertinente passa con **24 test case e 141407 asserzioni**, codice di uscita
**0**, sia per la compilazione sia per l'esecuzione. Comandi usati dalla
radice del repository, con il driver dedicato e i flag/macro/librerie della
build SU2 esistente:

```sh
python3 /private/tmp/su2-build-phase2-tests.py > /private/tmp/su2-phase2-extra-build.log 2>&1
/private/tmp/su2-phase2-tests '[primitive-geometry],[primitive-euler]' > /private/tmp/su2-phase2-extra-test.log 2>&1
```

Il driver usa Apple Clang 15 e C++17 e viene collegato alle librerie esistenti;
la build principale non è stata riconfigurata. Gli avvisi sono gli stessi
già riportati: pragma `-Wstringop-overflow` non riconosciuto da Clang,
`-lmpi`/`-lz` duplicati nella link line e bind TCP MPI negato dal sandbox
durante l'avvio. Quest'ultimo non impedisce l'esecuzione seriale, terminata
correttamente. Nessun test della selezione è stato saltato. L'intera suite
SU2, MPI distribuito e simulazioni del solver completo non sono stati
eseguiti perché esterni a questa verifica.

Le modifiche di questa verifica riguardano esclusivamente test e
documentazione: API e implementazione di `CAggregatedGeometry`, kernel e
solver conservativo rimangono invariati. Non sono stati creati commit o
push e non è stata avviata la FASE 3.
