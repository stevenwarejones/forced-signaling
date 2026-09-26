# Shot budgets under i.i.d. sampling

A: rigorous sufficient total N for at least 95% power; minimum integer allocation certified by this concentration bound, not the true minimum.
B: seeded Monte Carlo estimated 95%-power crossing (0.5% search resolution); independent validation power and pointwise exact 95% binomial interval.

| p | α | A scalar | A directional | B scalar: N; power [CI] | B directional: N; power [CI] |
|---|---|---:|---:|---|---|
| 0.90 | 0.05 | 14,134,528 | 14,134,528 | 6,294,272; 0.9492 [0.9460, 0.9522] | 6,045,808; 0.9500 [0.9469, 0.9530] |
| 0.90 | 0.01 | 15,386,496 | 15,386,496 | 7,182,352; 0.9568 [0.9539, 0.9596] | 6,911,888; 0.9573 [0.9544, 0.9601] |
| 0.90 | 0.001 | 17,087,024 | 17,087,024 | 8,343,264; 0.9516 [0.9485, 0.9545] | 8,042,896; 0.9528 [0.9498, 0.9557] |
| 0.92 | 0.05 | 3,763,072 | 3,763,072 | 1,675,728; 0.9478 [0.9446, 0.9508] | 1,609,568; 0.9541 [0.9511, 0.9570] |
| 0.92 | 0.01 | 4,096,384 | 4,096,384 | 1,904,160; 0.9533 [0.9503, 0.9562] | 1,824,160; 0.9472 [0.9440, 0.9502] |
| 0.92 | 0.001 | 4,549,120 | 4,549,120 | 2,221,232; 0.9550 [0.9520, 0.9578] | 2,132,400; 0.9492 [0.9460, 0.9522] |
| 0.95 | 0.05 | 1,263,120 | 1,263,120 | 560,016; 0.9512 [0.9481, 0.9541] | 537,808; 0.9547 [0.9517, 0.9575] |
| 0.95 | 0.01 | 1,375,008 | 1,375,008 | 633,776; 0.9488 [0.9457, 0.9518] | 609,616; 0.9470 [0.9438, 0.9501] |
| 0.95 | 0.001 | 1,526,976 | 1,526,976 | 739,616; 0.9510 [0.9480, 0.9540] | 715,760; 0.9594 [0.9566, 0.9621] |
| 0.98 | 0.05 | 625,872 | 625,872 | 276,256; 0.9522 [0.9492, 0.9551] | 265,248; 0.9595 [0.9566, 0.9621] |
| 0.98 | 0.01 | 681,296 | 681,296 | 312,688; 0.9490 [0.9459, 0.9520] | 300,720; 0.9490 [0.9459, 0.9520] |
| 0.98 | 0.001 | 756,592 | 756,592 | 364,992; 0.9507 [0.9476, 0.9537] | 351,696; 0.9513 [0.9482, 0.9542] |
| 1.00 | 0.05 | 436,528 | 436,528 | 191,824; 0.9507 [0.9477, 0.9537] | 184,144; 0.9564 [0.9535, 0.9592] |
| 1.00 | 0.01 | 475,184 | 475,184 | 218,096; 0.9526 [0.9496, 0.9555] | 209,744; 0.9580 [0.9551, 0.9607] |
| 1.00 | 0.001 | 527,712 | 527,712 | 254,560; 0.9571 [0.9543, 0.9599] | 244,256; 0.9516 [0.9485, 0.9545] |

The directional rejection region contains the scalar region for every dataset. Worst-case budgets coincide; Monte Carlo crossings have sampling uncertainty.

## Near the threshold

| p − p* | A total N (α=.01) |
|---|---:|
| 0.01 | 69,940,384 |
| 0.001 | 6,994,037,616 |
| 0.0001 | 699,403,760,432 |

p*=3−3√2/2. N grows as (p−p*)⁻²; no finite sufficient budget is certified at or below p*.

Seed 260926; 5,000 trials per search point; 20,000 independent validation trials per crossing.
Intervals are pointwise Monte Carlo uncertainty, not simultaneous guarantees or experimental significance bounds.
Full brackets, counts, environment: `data/shot_budget_results.json`.

## Optional nonuniform allocation

Numerically optimized feasible sufficient budgets for both tests. The condition is checked after rounding every setting count upward. Global optimality is not claimed; no nonuniform Monte Carlo study is claimed.

| p | α | Sufficient total N |
|---|---|---:|
| 0.90 | 0.05 | 12,889,291 |
| 0.90 | 0.01 | 13,978,259 |
| 0.90 | 0.001 | 15,458,485 |
| 0.92 | 0.05 | 3,431,560 |
| 0.92 | 0.01 | 3,721,478 |
| 0.92 | 0.001 | 4,115,566 |
| 0.95 | 0.05 | 1,151,859 |
| 0.95 | 0.01 | 1,249,171 |
| 0.95 | 0.001 | 1,381,449 |
| 0.98 | 0.05 | 570,744 |
| 0.98 | 0.01 | 618,963 |
| 0.98 | 0.001 | 684,505 |
| 1.00 | 0.05 | 398,085 |
| 1.00 | 0.01 | 431,719 |
| 1.00 | 0.001 | 477,430 |

All 16 positive integer setting counts are recorded in the JSON in lexicographic (x,y,z,w) order.
