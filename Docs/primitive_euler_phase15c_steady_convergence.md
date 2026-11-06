# FASE 15C — NACA 0012: campagna stazionaria con pseudo-tempo locale

## Esito e perimetro

**Convergenza stazionaria non raggiunta.** CFL=.05 mostra scenario B nel
tratto osservato; CFL=.1 e .2 mostrano scenario C, con risalite/oscillazioni
dopo la discesa iniziale. La campagna complessiva non valida il caso NACA.
Non sono sufficienti né il residuo inferiore alla FASE 15B né la positività.

Branch `feature/primitive-path-solver`, lavoro 15A/15B conservato nel working
tree. Nessun commit, push, sync o altra operazione Git di scrittura. Nessun
AOA=2°. Risultati del nativo 15A riutilizzabili ma **non ricalcolati**; nessun
confronto quantitativo Cp fra stati dichiarati stazionari.

Configurazione fisica invariata: mesh originale QuickStart, 5233 nodi e
10216 triangoli, Mach=.2, AOA=0°, gamma=1.4, rho=p=1, u=.2*sqrt(1.4), v=0;
NODE/AGGREGATED/NODE_POINT, WLS, MUSCL affine, limiter NONE, farfield
LOCAL_CHARACTERISTIC, TE debole Global_ID=199, REF_AREA=1. Tutte le prove
partono indipendentemente dal freestream, senza prosecuzioni/restart.
SHA256 mesh:
`9094b51c2628d3bb4c865d774b2308dbb59a213ad19c007cfcd7d57e7aaeebeb`.

Unica modifica produttiva di questa fase:
`SU2_CFD/src/output/CPrimitiveEulerOutput.cpp`, metodi
`SetHistoryOutputFields`/`LoadHistoryData` e
`SetVolumeOutputFields`/`LoadVolumeData`: aggiunte history della parete forte
e del farfield e componenti nodali del RHS nei VTU locali. Non modificati
operatore, stima spettrale, aggiornamento, BC, geometria, MUSCL o TE.
Nuovo orchestratore/postprocessor:
`UnitTests/run_primitive_euler_steady15c.py`.

## Budget, arresti e definizione delle norme

Budget configurabile di 480 s; campagna effettiva **431.94 s**, comprensiva
di preflight, scrittura e sorveglianza. Tre prove terminate normalmente al
limite di iterazioni: .05/3000, .1/3000, .2/6000. Non provato .4: le
oscillazioni ai CFL inferiori non giustificano un ulteriore aumento.

Sorveglianza ogni circa .5 s: arresto per budget, perdita di positività,
Mach nodale/ricostruito >=1, oppure quattro medie successive su 200 passi
ciascuna crescenti tutte di oltre il 5%. Nessuna prova ha attivato questi
arresti; le oscillazioni registrate non equivalgono a crescita monotona che
soddisfa quest'ultimo criterio. Un limite di iterazioni non è convergenza.

Si misura il **RHS spaziale corrente ricalcolato**, dopo la reazione dei
nodi forti, non la norma di tau_i*RHS:

    E_G = sqrt(sum_{i in G} V_i sum_{k=1}^4 RHS_ik^2 / sum_{i in G} V_i)
    e_G(n) = E_G(n)/E_G(0)
    d_G(n) = log10(E_G(0)/E_G(n))

Le regioni sono disgiunte: interno 4983 nodi, parete forte 199, TE 1,
farfield 50. Volumi rispettivi circa 1121.6812, .01154222, 9.73328e-8,
131.55777; volume totale 1253.2505. Ogni norma regionale usa il proprio
volume. Il TE ha anche norma assoluta, senza diluizione nel volume globale.
Il farfield inizialmente ha RHS zero: la sua normalizzazione iniziale è
lasciata vuota, non divisa per zero. Le componenti primitive sono quelle
adimensionali del caso; la somma delle quattro non è una norma di energia.

## Checkpoint e andamento aerodinamico

Valori assoluti; normalizzati e stato completo rho,u,v,p al TE in
`results15c/checkpoints.csv`. Tutte le prove hanno lo stesso stato iniziale:
Eglob=.0470270, Eint=.0260361, Ewall=13.19414, ETE=140.2453, CL=CD=0.

| CFL | Iter. | Eglob | Eint | Ewall | ETE | CL | CD |
|---:|---:|---:|---:|---:|---:|---:|---:|
| .05 | 300 | .0275565 | .0283484 | 2.08630 | 10.9326 | .00000305 | .683667 |
| .05 | 1000 | .0176420 | .0185143 | .694814 | 3.76278 | .00030591 | .398830 |
| .05 | 2000 | .0105513 | .0110523 | .465825 | 5.06873 | -.0146756 | .147923 |
| .05 | 3000 | .00744590 | .00778009 | .370533 | 4.52723 | -.0212537 | -.0542186 |
| .1 | 300 | .0225458 | .0236921 | .802077 | 4.50719 | -.00004458 | .508490 |
| .1 | 1000 | .0107988 | .0113149 | .468923 | 5.11289 | -.0146813 | .147883 |
| .1 | 2000 | .00552211 | .00568345 | .413118 | 9.73034 | -.00660229 | -.147724 |
| .1 | 3000 | .00687079 | .00710657 | .451027 | 24.2914 | .0194122 | -.0791987 |
| .2 | 300 | .0168016 | .0176572 | .593803 | 3.70689 | -.00044659 | .341817 |
| .2 | 1000 | .00578626 | .00594922 | .440975 | 9.36606 | -.00760663 | -.147800 |
| .2 | 2000 | .00541952 | .00563321 | .322012 | 13.1313 | .0165250 | .232007 |
| .2 | 3000 | .00931032 | .00952107 | .760951 | 52.5519 | .0122953 | -.0782238 |
| .2 | 5000 | .00610848 | .00626015 | .492165 | 9.39013 | .109436 | .129154 |
| .2 | 6000 | .00435403 | .00437464 | .439895 | 24.2570 | -.0774417 | -.223352 |

L'interno, inizialmente crescente a .05, diminuisce entro 1000 passi;
non rimane però in decadimento monotono ai CFL superiori. A .1 le medie
globali nelle finestre 2201–2400, 2401–2600, 2601–2800, 2801–3000 sono
.00461404, .00459677, .00526185, .00632209. A .2 le ultime cinque medie
su 200 passi sono .00555966, .00486002, .00609210, .00619343, .00495197:
un minimo finale favorevole non prova un decadimento asintotico.

| CFL | Decadi globali | Decadi interno | Decadi TE | Range CL ultimi 200 | Range CD ultimi 200 |
|---:|---:|---:|---:|---:|---:|
| .05 | .8004 | .5246 | 1.4911 | .001641 | .032934 |
| .1 | .8353 | .5639 | .7614 | .009880 | .030397 |
| .2 | 1.0335 | .7746 | .7621 | .051580 | .092483 |

CL si allontana da zero, CD cambia segno e nessuno dei due si stabilizza.
Questi sono valori di stati transitori in pseudo-tempo, non coefficienti
del profilo validati. La simmetria attesa non autorizza a imporre CL=0.

Stati TE finali:

| CFL | rho | u | v | p |
|---:|---:|---:|---:|---:|
| .05 | 1.00082645 | .20247573 | -.000691565 | 1.00119787 |
| .1 | 1.03638763 | .17236957 | .001760729 | 1.05133970 |
| .2 | 1.05028932 | .15512831 | -.001742011 | 1.07113716 |

## Costo, ammissibilità e passi locali

| CFL | Passi | Costo s | ms/passo | Decadi globali/s | min rho | min p | max M nodale | max M tracce |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| .05 | 3000 | 96.83 | 32.28 | .008266 | .967159 | .954327 | .238992 | .270820 |
| .1 | 3000 | 98.89 | 32.96 | .008447 | .966194 | .953088 | .242584 | .270820 |
| .2 | 6000 | 231.55 | 38.59 | .004463 | .908950 | .874896 | .250384 | .270820 |

Costi end-to-end con output e monitoraggio; non benchmark isolati del kernel.
Minimi/massimi coprono tutta la history, non soltanto il campo finale.
Zero tentativi rifiutati, nessuna riduzione di CFL; CFL effettivi .05/.1/.2.
Assenza di rifiuti non implica stabilità asintotica. Il CFL maggiore accelera
la prima discesa ma non migliora necessariamente la convergenza; le decadi/s
medie non descrivono il plateau e non vanno estrapolate come una legge.

Distribuzione finale tau_i (percentili nodali non pesati; valori VTU Float32):

| CFL | min = TE | p10 | p25 | mediana | p75 | p90 | max |
|---:|---:|---:|---:|---:|---:|---:|---:|
| .05 | 2.64069e-6 | 4.05081e-5 | 9.99959e-5 | 2.77340e-4 | .00141032 | .00805246 | .0241533 |
| .1 | 5.32190e-6 | 8.14332e-5 | 2.00506e-4 | 5.57924e-4 | .00283378 | .0160937 | .0483106 |
| .2 | 1.07014e-5 | 1.63212e-4 | 4.00084e-4 | .00111056 | .00565593 | .0322197 | .0966240 |

Non esiste un tempo fisico comune associato a questi passi locali.

## Localizzazione e limiti della diagnosi

All'inizio la parete forte pesa 72.50% nel quadrato della norma globale,
l'interno 27.43%, il TE .0691%. Alla fine:

| CFL | Interno | Parete forte | TE | Farfield | Efarfield assoluto |
|---:|---:|---:|---:|---:|---:|
| .05 | 97.716% | 2.281% | .00287% | .000165% | 2.94867e-5 |
| .1 | 95.750% | 3.969% | .0971% | .1846% | .000911087 |
| .2 | 90.351% | 9.401% | .2411% | .00678% | .000110632 |

Queste percentuali sono V_G*E_G²/(Vtot*Eglob²), non confronti fra norme
regionali normalizzate. **La mancata riduzione globale coinvolge soprattutto
l'interno**, per peso volumetrico; il TE mantiene un difetto locale elevato
e oscillante. Non è dimostrato che il TE generi le oscillazioni interne,
né che il farfield sia innocuo per la dinamica solo perché pesa poco nella
norma. I dati non distinguono definitivamente instabilità dell'iterazione
esplicita, risposta dell'operatore spaziale e accoppiamento parete/TE/interno.
Non sono stati introdotti Kutta, proiezioni aggiuntive o rimedi empirici.

La sola coda .05, 2000→3000, guadagna circa .151 decadi/1000 passi.
Mantenendo irrealisticamente quel tasso, arrivare a 3 decadi globali totali
richiederebbe circa altri 14600 passi (~470 s); per 4 circa altri 21200
(~685 s). È una stima condizionata, **non una previsione affidabile**:
TE e coefficienti sono già non stazionari, e gli altri CFL mostrano che il
decadimento può arrestarsi. Non autorizza una prosecuzione indefinita.

Il gate di postprocessing richiede almeno 3 decadi globali e interne,
ETE<1e-3 con riduzione relativa >=3 decadi, range CL/CD<1e-5 sugli ultimi
200 passi, uscita normale e snapshot finale presente. È un filtro prudente
di questa campagna, non una soglia universale di accuratezza. Nessun run lo
supera: `cp_comparison_status.txt` è False. Nessun Cp primitivo/nativo
quantitativo generato. Una successiva indagine deve prima separare cause
temporali e spaziali delle oscillazioni, usando i campi già salvati.

## Artefatti e verifiche riproducibili

Radice: `TestCases/primitive_euler_examples/naca0012_aoa0/results15c/`.

- `summary.csv`, `checkpoints.csv`, `campaign_execution.json` e
  `output_regression.json`: costi, norme, stati TE e controllo di invarianza.
- `cfl_0p05/`, `cfl_0p1/`, `cfl_0p2/`: configurazioni, log, history completa,
  `regional_checkpoints.csv`, `windows.csv`, `step_distribution.csv`,
  tabelle pressione/TE e restart diagnostici. Questi ultimi sono conservati,
  ma una prosecuzione locale da restart non è validata in questa fase.
- `iterations.pvd` in ciascuna cartella: checkpoint 0,300,1000,2000,3000;
  anche 5000,6000 per .2. Ascissa PVD = pseudo-iterazione. Altri VTU ogni
  300 passi e VTM ogni 1000 sono conservati. Campi rho, velocità, p, Mach,
  Global_ID, Dual_Volume, Local_Pseudo_Step, Pseudo_Iteration e quattro RHS.
  Nessun campo Physical_Time; eventuale TimeValue del framework non è tempo
  fisico. La sequenza punta direttamente ai VTU interni indicizzati,
  evitando blocchi VTM non indicizzati mescolati dal writer nativo.
- Grafici: `residuals.png`, `trailing_edge.png`, `regions.png`, `lift.png`,
  `drag.png`, tutti rispetto alla pseudo-iterazione; RHS globale e regioni
  non nulle normalizzati all'iniziale, TE anche assoluto.

Validati XML/PVD e riferimenti, payload binari VTU, numero nodi/elementi,
ID, indici di pseudo-iterazione, positività e Mach. Le cinque norme
ricostruite indipendentemente da RHS e volumi VTU coincidono con la history
entro rtol=3e-6 (precisione Float32); le quattro partizioni ricompongono il
quadrato globale. Non è stata effettuata una verifica GUI in ParaView.

Preflight di 50 passi a .05: restart iniziale/finale identico al caso smoke50
15B, massimo scarto **0.0**. La diagnostica aggiunta non cambia l'evoluzione.
Controllo CLI poligonale passato: 4 snapshot, 4 tabelle pressione, CL/CD
coerenti con le forze delle mezze facce. Build riuscita; nessuna ripetizione
dell'intera suite 15B. Il primo postprocessing ha incontrato soltanto un
errore di serializzazione JSON di un booleano NumPy, corretto convertendolo
in bool Python; ripetuto solo il postprocessing, **nessuna simulazione**.

```sh
./ninja -C build SU2_CFD/src/SU2_CFD
python3 UnitTests/run_primitive_euler_steady15c.py \
  --output-dir TestCases/primitive_euler_examples/naca0012_aoa0/results15c \
  --budget-seconds 480
# Per rigenerare tabelle/grafici SENZA simulazioni:
python3 UnitTests/run_primitive_euler_steady15c.py \
  --output-dir TestCases/primitive_euler_examples/naca0012_aoa0/results15c \
  --postprocess-only
python3 UnitTests/check_primitive_euler_polygon_cli.py \
  --executable build/SU2_CFD/src/SU2_CFD \
  --output-dir /private/tmp/primitive-phase15c-polygon
python3 -m py_compile UnitTests/run_primitive_euler_steady15c.py
git diff --check
```

L'orchestratore rifiuta una cartella di nuova campagna già esistente: per
ripetere i calcoli scegliere una nuova directory, senza sovrascrivere 15A/B/C.
