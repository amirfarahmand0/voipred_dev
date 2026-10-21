# Many model families on GUSTO-I. For each: nb_val() must give the same net
# benefit as dcurves::dca() at every threshold. Models that do not return
# probabilities by default have their own predict_risk() method; their raw,
# unconverted output (scores, classes) must be rejected.

check_model <- function(fit, val) {
  expect_same_as_dcurves(nb_val(val, fit, "day30", all_z), val)
}

# ---- 1. predict() already returns probabilities (default method) ----------

test_that("glm with probit and cloglog links agrees with dcurves", {
  s <- gusto_split()
  check_model(glm(gusto_formula, family = binomial(link = "probit"), data = s$dev), s$val)
  check_model(glm(gusto_formula, family = binomial(link = "cloglog"), data = s$dev), s$val)
})

test_that("Firth logistic regression (logistf) agrees with dcurves", {
  skip_if_not_installed("logistf")
  s <- gusto_split()
  check_model(logistf::logistf(gusto_formula, data = s$dev), s$val)
})

test_that("generalized additive model (mgcv::gam) agrees with dcurves", {
  skip_if_not_installed("mgcv")
  s <- gusto_split()
  fit <- mgcv::gam(day30 ~ s(age) + s(sysbp) + s(pulse) + Killip + pmi + miloc + sex,
                   family = binomial, data = s$dev)
  check_model(fit, s$val)
})

test_that("MARS (earth) with a binomial glm agrees with dcurves", {
  skip_if_not_installed("earth")
  s <- gusto_split()
  fit <- earth::earth(gusto_formula, data = s$dev, glm = list(family = binomial))
  check_model(fit, s$val)
})

test_that("classification tree (rpart) agrees with dcurves", {
  skip_if_not_installed("rpart")
  s <- gusto_split()
  fit <- rpart::rpart(factor(day30) ~ age + Killip + sysbp + pulse + pmi + miloc + sex,
                      data = s$dev, method = "class", cp = 0.002)
  check_model(fit, s$val)
})

test_that("conditional inference tree (partykit::ctree) agrees with dcurves", {
  skip_if_not_installed("partykit")
  s <- gusto_split()
  fit <- partykit::ctree(factor(day30) ~ age + Killip + sysbp + pulse + pmi + miloc + sex,
                         data = s$dev)
  check_model(fit, s$val)
})

test_that("gradient boosting (gbm) agrees with dcurves", {
  skip_if_not_installed("gbm")
  s <- gusto_split()
  set.seed(1)
  fit <- gbm::gbm(gusto_formula, data = s$dev, distribution = "bernoulli",
                  n.trees = 200, interaction.depth = 2)
  expect_same_as_dcurves(suppressMessages(nb_val(s$val, fit, "day30", all_z)), s$val)
})

test_that("k-nearest neighbours (kknn) agrees with dcurves", {
  skip_if_not_installed("kknn")
  s <- gusto_split()
  ff <- factor(day30) ~ age + Killip + sysbp + pulse + pmi + miloc + sex
  withr::local_package("kknn")   # kknn only predicts when attached
  fit <- kknn::train.kknn(ff, data = s$dev, kmax = 50)
  check_model(fit, s$val)
})

test_that("naive Bayes (naivebayes package) agrees with dcurves", {
  skip_if_not_installed("naivebayes")
  s <- gusto_split()
  fit <- naivebayes::naive_bayes(factor(day30) ~ age + Killip + sysbp + pulse + pmi + miloc + sex,
                                 data = s$dev)
  check_model(fit, s$val)
})

# ---- 2. models with their own predict_risk() method ----------------------

test_that("support vector machine (e1071::svm, probability = TRUE) agrees with dcurves", {
  skip_if_not_installed("e1071")
  s <- gusto_split()
  set.seed(1)
  fit <- e1071::svm(factor(day30) ~ age + Killip + sysbp + pulse + pmi + miloc + sex,
                    data = s$dev[1:2000, ], probability = TRUE)
  check_model(fit, s$val)
})

test_that("naive Bayes (e1071::naiveBayes) agrees with dcurves", {
  skip_if_not_installed("e1071")
  s <- gusto_split()
  fit <- e1071::naiveBayes(factor(day30) ~ age + Killip + sysbp + pulse + pmi + miloc + sex,
                           data = s$dev)
  check_model(fit, s$val)
})

test_that("neural network (nnet) agrees with dcurves, factor and 0/1 outcome", {
  skip_if_not_installed("nnet")
  s <- gusto_split()
  set.seed(1)
  fit_f <- nnet::nnet(factor(day30) ~ age + Killip + sysbp + pulse + pmi + miloc + sex,
                      data = s$dev, size = 3, decay = 0.1, maxit = 300, trace = FALSE)
  check_model(fit_f, s$val)
  set.seed(1)
  fit_n <- nnet::nnet(gusto_formula, data = s$dev, size = 3, decay = 0.1,
                      maxit = 300, entropy = TRUE, trace = FALSE)
  check_model(fit_n, s$val)
})

test_that("linear and quadratic discriminant analysis (MASS) agree with dcurves", {
  skip_if_not_installed("MASS")
  s <- gusto_split()
  ff <- factor(day30) ~ age + sysbp + pulse + Killip
  check_model(MASS::lda(ff, data = s$dev), s$val)
  check_model(MASS::qda(factor(day30) ~ age + sysbp + pulse, data = s$dev), s$val)
})

test_that("svm and naiveBayes ignore missing values in columns the model does not use", {
  skip_if_not_installed("e1071")
  s <- gusto_split()
  expect_true(any(!complete.cases(s$val)))     # NAs in unused GUSTO columns
  set.seed(1)
  fit <- e1071::svm(factor(day30) ~ age + Killip + sysbp + pulse, data = s$dev[1:2000, ],
                    probability = TRUE)
  expect_length(predict_risk(fit, s$val), nrow(s$val))
  bayes <- e1071::naiveBayes(factor(day30) ~ age + Killip + sysbp + pulse, data = s$dev)
  expect_length(predict_risk(bayes, s$val), nrow(s$val))
})

# ---- 3. models trained on a numeric matrix -------------------------------

test_that("glmnet (one lambda) and cv.glmnet agree with dcurves", {
  skip_if_not_installed("glmnet")
  m <- gusto_matrix_split()
  set.seed(1)
  check_model(glmnet::glmnet(m$x, m$y, family = "binomial", lambda = 0.02), m$val)
  check_model(glmnet::cv.glmnet(m$x, m$y, family = "binomial"), m$val)
})

test_that("LightGBM (objective binary) agrees with dcurves", {
  skip_if_not_installed("lightgbm")
  m <- gusto_matrix_split()
  fit <- lightgbm::lightgbm(data = m$x, label = m$y, nrounds = 100,
                            params = list(objective = "binary", learning_rate = 0.05,
                                          num_threads = 1, verbose = -1))
  check_model(fit, m$val)
})

test_that("XGBoost (binary:logistic) agrees with dcurves, both interfaces", {
  skip_if_not_installed("xgboost")
  m <- gusto_matrix_split()
  fit <- xgboost::xgb.train(params = list(objective = "binary:logistic", max_depth = 3,
                                          eta = 0.1, nthread = 1),
                            data = xgboost::xgb.DMatrix(m$x, label = m$y),
                            nrounds = 100, verbose = 0)
  check_model(fit, m$val)
  # xgboost() in version 3+ takes a data frame with factors, and a factor
  # outcome; the risk is P(day30 = 1) whatever the order of the levels
  skip_if(utils::packageVersion("xgboost") < "3.0.0")
  s <- gusto_split()
  vars <- c("age", "Killip", "sysbp", "pulse", "pmi", "miloc", "sex")
  fit_df <- xgboost::xgboost(s$dev[, vars], factor(s$dev$day30), nrounds = 100,
                             max_depth = 3, nthreads = 1)
  check_model(fit_df, s$val)
  fit_rev <- xgboost::xgboost(s$dev[, vars], factor(s$dev$day30, levels = c(1, 0)),
                              nrounds = 100, max_depth = 3, nthreads = 1)
  p_rev <- predict_risk(fit_rev, s$val)
  expect_equal(mean(p_rev), mean(s$val$day30), tolerance = 0.2)   # near prevalence, not 1 - it
  check_model(fit_rev, s$val)
})

test_that("CatBoost (Logloss) agrees with dcurves, numeric and categorical features", {
  skip_if_not_installed("catboost")
  m <- gusto_matrix_split()
  pars <- list(loss_function = "Logloss", iterations = 100, logging_level = "Silent",
               random_seed = 1, thread_count = 1, allow_writing_files = FALSE)
  fit <- catboost::catboost.train(catboost::catboost.load_pool(m$dev[, colnames(m$x)], label = m$y),
                                  params = pars)
  check_model(fit, m$val)
  # CatBoost also takes factors directly, without dummy coding
  s <- gusto_split()
  vars <- c("age", "Killip", "sysbp", "pulse", "pmi", "miloc", "sex")
  fit_cat <- catboost::catboost.train(catboost::catboost.load_pool(s$dev[, vars], label = s$dev$day30),
                                      params = pars)
  check_model(fit_cat, s$val)
})

# ---- 4. unconverted output is rejected -----------------------------------

test_that("svm without probability = TRUE, or its decision values, are rejected", {
  skip_if_not_installed("e1071")
  s <- gusto_split()
  ff <- factor(day30) ~ age + Killip + sysbp + pulse
  fit <- e1071::svm(ff, data = s$dev[1:2000, ])
  expect_error(nb_val(s$val, fit, "day30"), "probability = TRUE")
  vars <- c("age", "Killip", "sysbp", "pulse")
  dv <- attr(predict(fit, s$val[, vars], decision.values = TRUE), "decision.values")[, 1]
  expect_error(nb_val(s$val, fixed_model(dv), "day30"), "not probabilities")
})

test_that("log-odds from glmnet, gbm and probit glm are rejected as scores", {
  m <- gusto_matrix_split()
  fit_net <- glmnet::glmnet(m$x, m$y, family = "binomial", lambda = 0.02)
  link <- predict(fit_net, newx = as.matrix(m$val[, colnames(m$x)]))   # default: link
  expect_error(nb_val(m$val, fixed_model(as.vector(link)), "day30"), "not probabilities")
  skip_if_not_installed("gbm")
  s <- gusto_split()
  fit_gbm <- gbm::gbm(gusto_formula, data = s$dev, distribution = "bernoulli", n.trees = 100)
  link_gbm <- suppressMessages(predict(fit_gbm, s$val))                # default: link
  expect_error(nb_val(s$val, fixed_model(link_gbm), "day30"), "not probabilities")
})

test_that("XGBoost raw log-odds (binary:logitraw) are rejected", {
  skip_if_not_installed("xgboost")
  m <- gusto_matrix_split()
  fit <- xgboost::xgb.train(params = list(objective = "binary:logitraw", nthread = 1),
                            data = xgboost::xgb.DMatrix(m$x, label = m$y),
                            nrounds = 30, verbose = 0)
  expect_error(nb_val(m$val, fit, "day30"), "not probabilities")
})

test_that("glmnet with a whole lambda path is rejected", {
  m <- gusto_matrix_split()
  fit <- glmnet::glmnet(m$x, m$y, family = "binomial")
  expect_error(nb_val(m$val, fit, "day30"), "single lambda")
})

test_that("LightGBM raw scores and CatBoost raw values and classes are rejected", {
  m <- gusto_matrix_split()
  xv <- as.matrix(m$val[, colnames(m$x)])
  if (requireNamespace("lightgbm", quietly = TRUE)) {
    fit <- lightgbm::lightgbm(data = m$x, label = m$y, nrounds = 50,
                              params = list(objective = "binary", num_threads = 1, verbose = -1))
    raw <- predict(fit, xv, type = "raw")
    expect_error(nb_val(m$val, fixed_model(as.vector(raw)), "day30"), "not probabilities")
  }
  skip_if_not_installed("catboost")
  fit_cb <- catboost::catboost.train(catboost::catboost.load_pool(m$dev[, colnames(m$x)], label = m$y),
                                     params = list(loss_function = "Logloss", iterations = 50,
                                                   logging_level = "Silent", thread_count = 1,
                                                   allow_writing_files = FALSE))
  pool <- catboost::catboost.load_pool(m$val[, colnames(m$x)])
  raw_cb <- catboost::catboost.predict(fit_cb, pool, prediction_type = "RawFormulaVal")
  expect_error(nb_val(m$val, fixed_model(raw_cb), "day30"), "not probabilities")
  cls_cb <- catboost::catboost.predict(fit_cb, pool, prediction_type = "Class")
  expect_error(nb_val(m$val, fixed_model(cls_cb), "day30"), "look like predicted classes")
})

test_that("LDA discriminant scores are rejected", {
  skip_if_not_installed("MASS")
  s <- gusto_split()
  fit <- MASS::lda(factor(day30) ~ age + sysbp + pulse + Killip, data = s$dev)
  ld <- predict(fit, s$val)$x[, 1]
  expect_error(nb_val(s$val, fixed_model(ld), "day30"), "not probabilities")
})

test_that("predicted classes from rpart and naiveBayes are rejected", {
  skip_if_not_installed("rpart")
  skip_if_not_installed("e1071")
  s <- gusto_split()
  ff <- factor(day30) ~ age + Killip + sysbp + pulse + pmi + miloc + sex
  tree <- rpart::rpart(ff, data = s$dev, method = "class", cp = 0.002)
  cls <- as.numeric(as.character(predict(tree, s$val, type = "class")))
  expect_error(nb_val(s$val, fixed_model(cls), "day30"), "look like predicted classes")
  bayes <- e1071::naiveBayes(ff, data = s$dev)
  cls2 <- as.numeric(as.character(predict(bayes, s$val, type = "class")))
  expect_error(nb_val(s$val, fixed_model(cls2), "day30"), "look like predicted classes")
})

# ---- 5. matrix models need the same columns in the validation data -------

test_that("matrix models: missing or non-numeric feature columns give clear errors", {
  m <- gusto_matrix_split()
  fit <- glmnet::glmnet(m$x, m$y, family = "binomial", lambda = 0.02)
  v <- m$val[, setdiff(names(m$val), c("KillipII", "sexfemale"))]
  expect_error(nb_val(v, fit, "day30"), "missing 2 covariates.*KillipII, sexfemale")
  v2 <- m$val
  v2$age <- as.character(v2$age)
  expect_error(nb_val(v2, fit, "day30"), "not numeric: age")
})
