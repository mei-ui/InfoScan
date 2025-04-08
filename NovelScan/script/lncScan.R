#!/usr/bin/env Rscript
# lncScan.R - Predict lncRNA using SVM model

suppressPackageStartupMessages({
  library(caret)
  library(pROC)
  library(ROCR)
  library(dplyr)
  library(randomForest)
  library(e1071)
  library(xgboost)
  library(MLmetrics)
})

# 定义核心预测函数
lncScan <- function(feature, model,model_type,species) {
  message("\n+ Initiating prediction...")
  
  # 确保输入为数据框
  if (!is.data.frame(feature)) feature <- as.data.frame(feature)
  
  # 执行预测
  #prob.res <- predict(model, newdata = feature, type = "prob", decision.values = TRUE)
  if(model_type == "svm") {
  res <- predict(model, newdata = feature)
  
  # 构建结果数据框
  results <- data.frame(
    ID = rownames(feature),
    Prediction = res,
    #Coding.Potential = prob.res[, 1],
    stringsAsFactors = FALSE
  )
  }
  if(model_type == "XGB") {
  feature2 = as.matrix(feature)
  feature2 = feature2[,c("ORF.Max.Len"  , "ORF.Max.Cov" ,  "Fickett_Score", "conservation")]
  feature2 <- xgboost::xgb.DMatrix(data = feature2)
  res <- predict(model,feature2)
  if(species == "human") {
    pred_class <- ifelse(res > 0.42, 'non_coding', "coding")
  }
  if(species == "mouse") {
    pred_class <- ifelse(res > 0.31, 'non_coding', "coding")
  }
  #print(head(res))
  # 构建结果数据框
  results <- data.frame(
    ID = rownames(feature),
    Prediction = pred_class,
    #Coding.Potential = prob.res[, 1],
    stringsAsFactors = FALSE
  )
  }
  message("+ Prediction completed successfully.")
  return(results)
}

# 主函数处理参数和流程
main <- function() {
  # 解析命令行参数
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 5) {
    stop("Usage: Rscript lncScan.R <species> <input_file> <output_file> <workdir> <model_type>\n",
         "  species:     'human' or 'mouse'\n",
         "  input_file:  input features file (tab-separated)\n",
         "  output_file: output results file\n",
         "  workdir:  input workdir\n",
         "  model_type:  input model type")
  }
  
  species <- tolower(args[1])
  input_file <- args[2]
  output_file <- args[3]
  workdir <- args[4]
  model_type <- args[5]
  
  # 验证物种参数
  if (!species %in% c("human", "mouse")) {
    stop("Invalid species. Valid options: 'human' or 'mouse'")
  }
  
  # 加载模型文件
  model_file <- paste0(workdir,"/",species, "_",model_type,".mod")
  if (!file.exists(model_file)) {
    stop("Model file not found: ", model_file)
  }
  message("\n+ Loading ", species, " model (", model_file, ")...")
  if(model_type == "XGB") {
    model <- xgboost::xgb.load(model_file)
  }
  if(model_type == "svm") {
    model <- readRDS(model_file)
  }
  
  # 读取输入文件
  message("\n+ Reading input features from: ", input_file)
  feature <- tryCatch(
    {
      df <- read.delim(input_file, stringsAsFactors = FALSE, check.names = FALSE)
      df[is.na(df)] <- 0  # 替换NA值为0
      rownames(df) <- df$gene
      colnames(df) <- gsub("_", "\\.", colnames(df))
      colnames(df)[colnames(df) == "Fickett.Score"] <- "Fickett_Score" 
      df <- subset(df , select = -c(gene))
      df
    },
    error = function(e) stop("Error reading input file: ", e$message)
  )
  
  # 执行预测
  results <- lncScan(feature, model,model_type,species)
  
  # 写入输出文件
  message("\n+ Writing results to: ", output_file)
  write.table(results, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message("\n+ Analysis complete. Results saved to: ", normalizePath(output_file))
}

# 执行主程序
if (sys.nframe() == 0) {
  main()
}
