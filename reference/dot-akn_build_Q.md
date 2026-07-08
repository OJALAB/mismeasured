# Build inverse Q matrix and baseline probability vector

Q = \[p_1 - p_0, ..., p\_(K-1) - p_0\] is (K-1) x (K-1) where p_l =
(Pi\[2, l+1\], ..., Pi\[K, l+1\])^T = E\[u \| Z = l\] and u is the
(K-1)-dim dummy encoding of Z_hat with baseline 0.

## Usage

``` r
.akn_build_Q(Pi, K)
```
