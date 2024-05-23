install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
BiocManager::install("VariantAnnotation", force = TRUE)
Sys.setenv(OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJjaGVuc2g5NTA5MjVAZm94bWFpbC5jb20iLCJpYXQiOjE3MTQ1NjMyMjgsImV4cCI6MTcxNTc3MjgyOH0.IVkiV_Q1WK8L8oOfHLZCVxkG7lHKD6niNQ-RPcyvEaGlo7TgnoJ4TM50ZeSLOGB1qHIPvStM1xxoqyYHjiQnNR4F_wyNHtmEZpl_D9zrsN91BUexpyK7IcnKG0dW-xSQtXIluWPkU4hMCzrb0R6FxF-qhuTfnM0ahTDwoMUZghC2k7IGWQCDknF3uANc7Dy_sTlU526nsvlAxxN9AFovK53HR_PRz66_qAGZLlmh4sG0B2zLGmGBBmyO0mUsBoUHCcNA6KixV9laWWNIXjVG5lUsVXwiM75ASYScKVb_iyW4EHx6dwRVx8aiB6ROWW5mZjn5OZdA7D_9nKshMa-3Ng")

library(dplyr)
library(ieugwasr)
Sys.setenv(OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJjaGVuc2g5NTA5MjVAZm94bWFpbC5jb20iLCJpYXQiOjE3MTYyNTUwMjUsImV4cCI6MTcxNzQ2NDYyNX0.RuclNKl9lcBRPfKuTdBGvlIw3r6jf5M1qsorbpLDjvb6KxKIS8tt6LJGxXz43TQfO9U_mXU7dKOVe7hVzCbmVJiR8E9pXV94huTWsaOpsQnoqZ0WEXY60D7RejPaXNsIHoO0AY7yGkTYxLD2CZAGSlA67ZJMOgVTTqANSVNNtydZEhIP5gylvMbJtHI4i4-rwVjp-WfQSRlU0voQli_b_QG4VcvBKN8-OFppYifL-slSOobzXwsRTzP2YCEybxnTyxqNbOWC2rGWZjJq3L7_TKb9D9MTUr-fJ-dVkCZvJNtIY7rEZ4kHgx5XNS11gvDy_nWqYANNcpRwmOF0dVMZug")
# 查看当前token
ieugwasr::get_opengwas_jwt()
# 返回用户名、账号类型，token过期时间
ieugwasr::user()

library(TwoSampleMR)
library(VariantAnnotation)
library(gwasglue2)
ebi-a-GCST90018803
ukb-b-14960
exposure_dat<-extract_instruments(outcomes = 'ukb-b-13803',p1 = 5e-08)

exposure_dat<- extract_instruments(
  outcomes='ukb-b-6235',
  clump=TRUE, r2=0.001,
  kb=10000,access_token= NULL
)
write.csv(exposure_dat,file = "exposure_dat（胆囊切除术）.csv")
exposure_dat<-read.csv('exposure_dat（胆囊切除术去混杂）.csv',header = T)

MR_F <- function(sample_size,num_IVs,r_square){
  numberator <- r_square * (sample_size - 1 - num_IVs)
  denominator <- (1 - r_square) * num_IVs
  f <- numberator / denominator
  return(f)
}
my_PVE <- function(beta, se_beta, maf, N){
  numberator <- 2 * beta^2 * maf * (1 - maf)
  denominator <- 2 * beta^2 * maf * (1 - maf) + se_beta^2 * 2 * N * maf * (1 - maf)
  return(numberator/denominator)
}

exposure_dat$pve <- my_PVE(beta = exposure_dat$beta.exposure, 
                           se_beta = exposure_dat$se.exposure,
                           maf = exposure_dat$eaf.exposure,
                           N = exposure_dat$samplesize.exposure)


bmi<-system.file("exposure.csv",package="TwoSampleMR")  
bmi_exp_dat<-read_exposure_data(filename = bmi,sep = ",",snp_col ="SNP",beta_col="beta",se_col="standard_error",effect_allele_col="effect_allele",other_allele_col="other_allele",clump = TRUE)

outcome_dat<-extract_outcome_data(snps=exposure_dat$SNP,outcomes="ebi-a-GCST90018583")

dat <- harmonise_data(exposure_dat, outcome_dat)

mrResult=mr(dat)

mrTab=generate_odds_ratios(mrResult)

heterTab=mr_heterogeneity(dat)

pleioTab=mr_pleiotropy_test(dat)

mr_scatter_plot(mrResult, dat)


res_single=mr_singlesnp(dat)      

mr_forest_plot(res_single)

mr_funnel_plot(singlesnp_results = res_single)

mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
