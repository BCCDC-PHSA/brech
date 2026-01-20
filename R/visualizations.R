#' Plot a Contact Matrix
#'
#' This function visualizes a contact matrix using a heatmap.
#'
#' @param matrix_data A square matrix where rows and columns represent different groups, and values indicate contact rates.
#' @param levels list of labels to order matrix
#' @return A ggplot2 object representing the heatmap of the contact matrix.
#'
#' @examples
#' contact_matrix <- matrix(runif(25, 0, 1), nrow = 5, dimnames = list(LETTERS[1:5], LETTERS[1:5]))
#' plot_contact_matrix(contact_matrix)
#'
#' @export
plot_contact_matrix <- function(matrix_data,levels=NULL){
  x <- y <- fill <- value <- NULL
  g <- dplyr::as_tibble(matrix_data,rownames="x") |>
    tidyr::pivot_longer(-x, names_to = "y", values_to = "value")
  if(!is.null(levels)){
    g <- g |>
      dplyr::mutate(
        x = forcats::fct_relevel(x,levels),
        y = forcats::fct_relevel(y,levels)
      )
  }

  g <- g |>
    ggplot2::ggplot(ggplot2::aes(x, y, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "white", high = "blue") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(fill="Contact rate",x="",y="")
  return(g)
}

#' Add Labeled Bars to a ggplot Object by Age Group
#'
#' This function enhances a ggplot object by adding bars with labels indicating
#' positive values, formatted with a steel blue fill color. It also adjusts the
#' y-axis scaling and applies a classic theme.
#'
#' @param g A ggplot object to which the bars and labels will be added.
#' @param label the column name of the y variable
#'
#' @return A modified ggplot object with labeled bars grouped by age.
#'
#' @examples
#' library(ggplot2)
#' data <- data.frame(
#'   age_group = c("0-18", "19-35", "36-50", "51+"),
#'   positive = c(10, 25, 30, 15)
#' )
#' g <- ggplot(data, aes(x = age_group, y = positive))
#' g <- add_labeled_bars_age_group(g)
#' print(g)
#'
#' @export
add_labeled_bars_age_group <- function(g,label="positive"){
  .data <- NULL
  g <- g +
  ggplot2::geom_bar(stat="identity",fill="steelblue") +
    ggplot2::geom_text(colour = "white", size = 3.5,
              ggplot2::aes(label = .data[[label]]),
              position=ggplot2::position_stack(vjust=0.5)) +
    ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::theme_classic() +
    ggplot2::labs(x="age group")

  return(g)
}
