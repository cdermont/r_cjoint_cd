# cjoint_cd
some custom functions based on r::cjoint

Changes: 

`plot.amce2` based on the original `plot.amce`, with the small change that `print(p)` is the last argument (similar to older versions), which allows to save the plot as an R object and reuse the data within. Combine results from multiple models, control colors and facets more explicitly. You'll find the data in `object$plot$data`.

`use.qualtrics` introduced based on `read.qualtrics`, which takes an R object as argument to reshape data into cjoint long format. Recode your data before transforming it to cjoint format, or use the same function and setup to transform data not from Qualtrics.

`amce.lmer` which applies `lme4::lmer` instead of clustered OLS. LMER results are stored in the object with the results unter `object$lmer.prof` and `object$lmer.full`, everything else is similar to amce. Get random effects per individual instead of clustered errors through `data.frame(ranef(object$lmer.prof)$respondentid)`, for example. Define hierarchy as in lmer, e.g., `hierarchy = "(1 | respondentid)"`.

