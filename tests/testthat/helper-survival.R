# test-km_ratio.R calls survival::survfit()/Surv() unqualified; the package
# itself calls them namespaced (survival:: prefix) so this isn't needed by
# the package code, only by the tests.
library(survival)
