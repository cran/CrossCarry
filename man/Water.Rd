% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Water.R
\docType{data}
\name{Water}
\alias{Water}
\title{Water student data of crossover design}
\format{
A data frame with 214 rows and 10 columns:
\describe{
\item{ID}{The unique identifier of each of the students}
\item{Age}{The age in years of each of the students}
\item{LCS}{Letter Cancel incorrect (raw score)}
\item{LCC}{Letter Cancel correct (raw score)}
\item{LCI}{Letter Cancel score (LCC-LCI)}
\item{sex}{ f=female, 1=male}
\item{school}{school indicator A or B}
\item{Treatment}{Condition indicator 0=Control 1=Water}
\item{Period}{date of visit}
\item{Carry_Agua}{Carry indicator 0=First period, 1=Water in the first period,
-1 = Control in the first period}
}
}
\source{
<https://journals.plos.org/plosone/article/authors?id=10.1371/journal.pone.0210568>
}
\usage{
Water
}
\description{
A pilot study to investigate the impact of providing supplementary drinking
water on the cognitive performance of pupils of two school grades (5 and 6)
in water-scarce schools in rural Mali.
47 students were assigned to take the control treatment (normal conditions)
on the first day and receive the treatment (controlled hydration)
on the second day. 60 received the treatments in reverse
(Hydration the first day and control the second day).
One part of this test assesses visual attention.
This test assesses visual attention. Pupils were given a grid containing
target letters randomly dispersed among non-target letters and were given
one minute to cross out as many target letters as possible.
Scores were calculated by subtracting the number of non-target letters
identified from the number of target letters identified; the maximum
test score was 38.
}
\keyword{datasets}
