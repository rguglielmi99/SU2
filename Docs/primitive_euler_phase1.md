# Euler 2D primitivo — FASE 1

Riferimento: specifica v1.1 del 20 ottobre 2026. Per il primo solver, come
concordato, si procede con Eulero esplicito, AGGREGATED e NODE. SSP-RK2,
SEGMENTED e DUAL_MEAN sono rinviati. Eulero esplicito è del primo ordine nel
tempo: la dimostrazione dell'ordine due complessivo resta una fase successiva.

## Kernel e convenzioni

`SU2_CFD/include/numerics/flow/primitive_euler_2d.hpp` contiene funzioni inline,
indipendenti da mesh, configurazione, solver e MPI. Il kernel iniziale usa
`double`; non è ancora predisposto o verificato per i tipi AD di SU2.

- `State`: quattro componenti nell'ordine **rho, u, v, p**. Non passare il
  vettore `Primitive` di `CEulerVariable`, che ha un layout diverso.
- `PhysicalMatrix(W,n,gamma)`: matrice `nx*Ax(W)+ny*Ay(W)`.
- `PathInverseDensity(rhoL,rhoR)`: integrale esatto di `1/rho` lungo il path
  lineare. Usa `log1p` vicino a densità uguali e restituisce `1/rho` per
  densità identiche. Gli estremi sono ordinati internamente per rendere
  simmetrica anche la valutazione in virgola mobile.
- `PathMatrix(L,R,n,gamma)`: integrale esatto della matrice lungo il path.
  Usa medie aritmetiche di rho, u, v e p e la media logaritmica di `1/rho`.
- `PathSpectrum`: autovalori `un-c|n|, un, un, un+c|n|`, dove
  `c=sqrt(gamma*p_mean*beta)`. Non confondere questa velocità di interfaccia
  con il suono fisico `sqrt(gamma*p/rho)` in un singolo stato.
- `AbsolutePathMatrix`: valore assoluto spettrale tramite proiettori.
  Con `K=(A-un*I)/(c|n|)`, i proiettori sono `I-K²`, `(K²-K)/2`, `(K²+K)/2`.
  L'autovalore doppio viene quindi trattato senza inversione di autovettori.
- `ComputeFluctuations(L,R,S,gamma)`: riceve la **normale integrata** di SU2,
  la normalizza e applica `|S|` una sola volta. Restituisce `minus` e `plus`.

Il futuro assemblaggio dovrà eseguire:

```cpp
const auto d = PrimitiveEuler2D::ComputeFluctuations(Wi, Wj, Sij, gamma);
// Aggiungere d.minus a Ri e d.plus a Rj, entrambi con segno positivo.
// Aggiungere separatamente volume e bordo, poi applicare -R/Vi.
```

Una faccia di area nulla restituisce fluttuazioni nulle. Densità o pressione
non positive, gamma <= 1 e dati non finiti sono rifiutati con
`std::invalid_argument`. Non sono applicati entropy fix, clipping o fallback.
Questi controlli non sostituiscono la positività della ricostruzione e degli
aggiornamenti temporali. Il kernel assume scale numeriche rappresentabili in
double; non offre garanzie contro overflow/underflow per stati estremi.

## Test e riproduzione

`UnitTests/SU2_CFD/numerics/primitive_euler_2d_tests.cpp` è registrato nel
driver Catch2 di SU2, sotto il tag `[primitive-euler]`. I test coprono:

- matrici contro le quattro equazioni primitive espanse;
- beta, limite a densità uguali e simmetria;
- matrice del path contro quadratura Simpson indipendente;
- azione della matrice e del suo valore assoluto sui quattro campi
  caratteristici, incluse le due direzioni dell'autovalore doppio;
- consistenza `D- + D+ = A_path * DeltaW * area`;
- inversione di stati/normale e scala dell'area;
- stato costante, area nulla, propagazione supersonica e stati sonici;
- limite della matrice a stati uguali, contatto stazionario e rotazioni;
- rifiuto degli input non ammissibili.

Per eseguire solo questi test senza ricompilare SU2, dalla radice del repository:

```sh
cat > /private/tmp/su2-primitive-test-main.cpp <<'EOF'
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
EOF
c++ -std=c++11 -O1 -Wall -Wextra -Wpedantic -I externals/catch2 \
  /private/tmp/su2-primitive-test-main.cpp \
  UnitTests/SU2_CFD/numerics/primitive_euler_2d_tests.cpp \
  -o /private/tmp/su2-primitive-tests
/private/tmp/su2-primitive-tests '[primitive-euler]'
```

Con un driver SU2 già compilato con `-Denable-tests=true`, si può invece usare:

```sh
build/UnitTests/test_driver '[primitive-euler]'
```

## Lavoro successivo

Questo kernel non registra ancora un solver eseguibile in SU2. Restano da
implementare geometria ausiliaria necessaria, classi solver/variabili,
configurazione, gradienti e limiter, residuale di volume, bordi, output e
avanzamento temporale.

NODE non elimina la distinzione fra valore nodale e media duale. Per mantenere
MEAN_COEFF del PDF occorre ancora il centroide duale, anche senza DUAL_MEAN:
`Wmean = Wi + Glim * (xc_i - xi)`. Usare direttamente `Wi` nei coefficienti
del volume sarebbe la variante NODE_POINT, da dichiarare separatamente.
La beta di interfaccia non sostituisce la quadratura del termine di volume.
