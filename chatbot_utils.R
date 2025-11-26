# AI Chatbot Integration Module
# Uses Azure OpenAI / HuggingFace LLM for interactive analysis

library(httr2)
library(jsonlite)
library(tidyverse)

#' Initialize LLM Connection
#' Supports Azure OpenAI and HuggingFace endpoints
initialize_llm <- function(provider = "azure", endpoint = NULL, api_key = NULL) {
  message(sprintf("Initializing %s LLM...", provider))
  
  config <- list(
    provider = provider,
    endpoint = endpoint %||% Sys.getenv(sprintf("%s_ENDPOINT", toupper(provider))),
    api_key = api_key %||% Sys.getenv(sprintf("%s_API_KEY", toupper(provider))),
    model = Sys.getenv(sprintf("%s_MODEL", toupper(provider)), 
                       if (provider == "azure") "gpt-4" else "meta-llama/Llama-2-7b-hf")
  )
  
  return(config)
}

#' Send Query to LLM with Context
#' @param query User question
#' @param context Biological context (cell types, markers, DE genes)
#' @param config LLM configuration
#' @param llm_config List with provider, endpoint, api_key, model
query_llm <- function(query, context = "", llm_config) {
  
  system_prompt <- "You are an expert bioinformatician specializing in single-cell RNA-seq analysis. 
  Provide accurate, scientifically grounded responses about cell types, gene expression, 
  and immune cell biology. Always cite relevant biological knowledge."
  
  messages <- list(
    list(role = "system", content = system_prompt),
    list(role = "user", content = sprintf("%s\n\nContext: %s", query, context))
  )
  
  if (llm_config$provider == "azure") {
    response <- call_azure_openai(messages, llm_config)
  } else if (llm_config$provider == "huggingface") {
    response <- call_huggingface_api(messages, llm_config)
  } else {
    stop("Unknown LLM provider")
  }
  
  return(response)
}

#' Call Azure OpenAI API
call_azure_openai <- function(messages, config) {
  
  request_body <- list(
    messages = messages,
    temperature = 0.7,
    max_tokens = 500,
    top_p = 0.95
  )
  
  response <- request(paste0(config$endpoint, "/chat/completions")) %>%
    req_headers("api-key" = config$api_key) %>%
    req_method("POST") %>%
    req_body_json(request_body) %>%
    req_perform()
  
  if (resp_status(response) == 200) {
    content <- resp_body_json(response)
    return(list(
      response = content$choices[[1]]$message$content,
      tokens_used = content$usage$total_tokens,
      status = "success"
    ))
  } else {
    return(list(response = "Error calling API", status = "error"))
  }
}

#' Call HuggingFace Inference API
call_huggingface_api <- function(messages, config) {
  
  # Format messages as conversation
  prompt <- paste(
    sapply(messages, function(m) sprintf("%s: %s", m$role, m$content)),
    collapse = "\n"
  )
  
  request_body <- list(
    inputs = prompt,
    parameters = list(
      max_new_tokens = 500,
      temperature = 0.7,
      top_p = 0.95
    )
  )
  
  response <- request("https://api-inference.huggingface.co/models/") %>%
    req_headers(Authorization = sprintf("Bearer %s", config$api_key)) %>%
    req_method("POST") %>%
    req_body_json(request_body) %>%
    req_perform()
  
  if (resp_status(response) == 200) {
    content <- resp_body_json(response)
    return(list(
      response = content[[1]]$generated_text,
      status = "success"
    ))
  } else {
    return(list(response = "Error calling API", status = "error"))
  }
}

#' Create RAG Context from Gene Expression Data
#' Extracts relevant information for LLM context
create_rag_context <- function(seurat_obj, query_genes, top_n = 5) {
  
  # Get differentially expressed genes
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.1)
  
  # Get gene expression in relevant clusters
  context_data <- list(
    top_genes = markers %>% 
      group_by(cluster) %>%
      slice_head(n = top_n) %>%
      pull(gene) %>%
      unique(),
    cell_types = unique(Idents(seurat_obj)),
    cluster_sizes = table(Idents(seurat_obj))
  )
  
  # Format as text for LLM
  context_text <- sprintf(
    "Dataset: %d cells, %d genes\nCell types: %s\nTop marker genes: %s",
    ncol(seurat_obj),
    nrow(seurat_obj),
    paste(context_data$cell_types, collapse = ", "),
    paste(context_data$top_genes[1:min(10, length(context_data$top_genes))], collapse = ", ")
  )
  
  return(context_text)
}

#' Generate Cell Type Predictions using LLM
#' @param markers Marker genes per cluster
#' @param seurat_obj Seurat object
predict_celltypes_llm <- function(markers, llm_config) {
  
  # Format marker genes for LLM
  markers_text <- markers %>%
    group_by(cluster) %>%
    slice_head(n = 5) %>%
    summarise(genes = paste(gene, collapse = ", "), .groups = "drop") %>%
    mutate(prompt = sprintf("Cluster %s: %s", cluster, genes)) %>%
    pull(prompt) %>%
    paste(collapse = "\n")
  
  query <- sprintf(
    "Based on these marker genes, predict the cell type for each cluster:\n%s\nProvide cell type annotations.",
    markers_text
  )
  
  response <- query_llm(query, "", llm_config)
  return(response)
}

#' Store conversation history for context
chat_memory <- new.env()
chat_memory$history <- list()

add_to_memory <- function(role, content) {
  chat_memory$history[[length(chat_memory$history) + 1]] <- 
    list(role = role, content = content, timestamp = Sys.time())
}

get_memory <- function(n_messages = 10) {
  messages <- tail(chat_memory$history, n_messages)
  return(lapply(messages, function(m) list(role = m$role, content = m$content)))
}
