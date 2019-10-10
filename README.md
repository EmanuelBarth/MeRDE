# MeRDE
This is MeRDE, a statisticak tool to analyze microRNA expression data received from RNA-Seq libraries.

## MeRDE: A new statistical model to infer differential expression of small RNAs from read counts
First, we assume that for every sample \(y\) the amount of reads assigned to a any gene \(x\) can be modeled by a gamma distributed random variable \(G_{xy}\):

\[G_{xy}\sim\Gamma(P_{xy}, B_{xy}) ,\]

with a shape parameter \(P_{xy} > 0\) and a scale paramter \(B_{xy} > 0\), resulting in non-negative read count values. Thus, the respective probability density function \(f_{xy}\) is:

\[f_{xy}(z) = \begin{cases}
    \dfrac{B_{xy}^{P_{xy}}}{\Gamma(P_{xy})} \cdot z^{P_{xy}-1} \cdot e^{-B_{xy}z}, & z > 0 \\
    0, & z \le 0
    \end{cases},\]

where \(\Gamma(P_{xy})\) is the gamma function, evaluated at \(P_{xy}\). Figure [fig.merde:gamma<sub>p</sub>df] displays examples of probability density plots of gamma distributions.

![image](https://github.com/EmanuelBarth/MeRDE/gamma_pdf.pdf) [fig.merde:gamma<sub>p</sub>df]

The mean \(\mu_{xy}\) and variance \(\sigma_{xy}\) of the random variable \(G_{xy}\) are defined as:

\[\mu_{xy} = P_{xy} \cdot B_{xy}\] \[\sigma_{xy} = P_{xy} \cdot B_{xy}^{2}\]

Based on this, we further assume that the estimated expression strength \(E_{xy}\) of gene \(x\) in sample \(y\) is the product of our modeled gene count and a library-specific size factor \(S_y\):

\[E_{xy} = G_{xy} \cdot S_y\]

When fitting our model to data, we expect this data to be a \(x \times y\) count matrix, where \(x = 1 \dots n\) is the number of genes, \(y = 1 \dots m\) the number of samples and \(g_{xy}\) the count value of gene \(x\) in sample \(y\). Library size factors are necessary to normalize RNA-Seq libraries in respect to their sequencing depths, or in other words, to equalize the total amount of reads between sequencing libraries that are compared. The calculation of the size factor estimators \(\hat{s}_{y}\) is based on all libraries \(y = 1 \dots m\) that are part of a comparison and is performed similarly as suggested by :

\[\hat{s}_{y} = \underset{x}{median} \dfrac{g_{xy}}{(\prod_{i = 1}^{m} g_{xi})^{1/m}}\].

Subsequently, we can estimate the average expression strength \(e_{x,c(A)}\) of any gene \(x\) of some condition \(A\) with \(c(A)\) being all sampled read libraries of condition \(A\) and \(|c(A)| = r \le m\):

\[e_{x,c(A)} = \dfrac{\sum_{y \in c(A)} e_{xy}}{r} = \dfrac{\sum_{y \in c(A)} g_{xy} \cdot \hat{s}_{y}}{r}.\]

Since we do not know the shape and rate paramters \(P_{xy}\) and \(B_{xy}\) of our gene expression strengths \(E_{xy}\), we have to estimate them based on the given data. Fortunately, there exist good maximum-likelihood estimators for both parameters (within 1.5% of the correct value) \(\hat{P}_{xy}\) and \(\hat{B}_{xy}\), that can be calculated efficiently :

\[\hat{P}_{xy} \approx \dfrac{3 - t + \sqrt{(t - 3)^2 + 24t}}{12t}\]

with

\[t = ln\begin{pmatrix}\dfrac{\sum_{y \in c(A)}g_{xy}}{r}\end{pmatrix} - \begin{pmatrix}\dfrac{\sum_{y \in c(A)}ln(g_{xy})}{r}\end{pmatrix}\]

and

\[\hat{B}_{xy} = \dfrac{1}{r\hat{P}_{xy}} \cdot  \sum_{y \in c(A)}g_{xy}.\]

### Building gene expression clusters

With the above introduced random variables \(E\) and \(G\) as well as both parameter estimators \(\hat{P}\) and \(\hat{B}\) it is possible to model the gene expression strength of any given miRNA gene from small RNA-Seq data. However, the performances of the shape and scale estimators strongly depend on the available amount of data, *i.e.*, the number of independently measured count values for every gene in each condition. The higher the amount of biological replicates, the better perform \(\hat{P}\) and \(\hat{B}\). As previously mentioned in Section 4.4.1, the number of biological replicates per investigated condition for RNA-Seq experiments is typically limited, seldom exceeding five or more sequenced samples. To overcome this obstacle, we assume that genes sharing a similar mean expression strength also share a similar dispersion. With this assumption it is possible to cluster genes based on their mean expression strengths and use the available count data of all genes within each cluster respectively to estimate their common shape and scale parameters \(\hat{P}\) and \(\hat{B}\) (see Figure [fig.merde:cluster]).

![image](./figs/merde/clusters1) [fig.merde:cluster]

For every gene \(x\) we calculate its expression strength mean \(\mu_{x}\) as the empirical mean read count of all replicates of gene \(x\). And we define the *expression range* of every gene \(x\) as the interval \(R_x\):

\[R_x = [\mu_{x} - \alpha\sigma_x, \mu_{x} + \alpha\sigma_x] ,\]

with \(\sigma_x\) being the empirical standard deviation of the gene expression of gene \(x\) and \(\alpha\) an additional scaling factor. At the moment, \(\alpha\) can be assumed to be a strictly positive value and its precise function will be discussed later. Again, for every gene \(x\) we can now determine its expression strength cluster \(C(x)\), which includes all genes \(w\) whose expression mean \(\mu_y\) lies within \(R_x\). This approach is performed for each condition and their respective replicates separately, allowing to estimate \(P_x\) and \(B_x\) on the aggregated counts in \(C(x)\).

### Outlier detection and expression cluster correction

Outliers (*i.e.*, singular extreme values often deriving from experimental errors) within count data represent sources of potentially great bias and have to be eliminated. But they are hard to identify if only few replicates are available, thus, the aggregation of count values for each gene based on the introduced expression clusters can also help in that matter. However, because expression clusters are calculated using mean and standard deviation of gene counts, they can be strongly affected by outliers themselves. In the presence of one or more outlier values in a given gene \(x\) its associated expression range \(R_x\) is greater than it should be, leading to an inclusion of genes \(w\) into the expression cluster \(C(x)\) that do not fit to the true mean \(\mu_x\). Accordingly, the subsequent estimation of the shape and scale parameters \(P_x\) and \(B_x\) can be strongly biased towards the outlier values. We can still use the information provided by the expression clusters to detect the presence of outliers within the examined dataset. For that, we developed two criteria that mark gene clusters which are potentially affected by outliers:

1.  Single count values of gene \(x\) are outside the expression interval spanned by its expression cluster \(C(x)\).

2.  The number of genes \(w\) assigned to expression cluster \(C(x)\) is higher than would be expected considering its mean value \(\mu_{C(x)}\).

Similarly to the gene expression range \(R_x\), the interval \(R_{C(x)}\) spanned by an expression cluster is defined by its mean and standard deviation, but without utilizing an additional scaling parameter \(\alpha\).
If only single outliers are present within an expression cluster \(C(x)\) its dispersion is not altered strongly, making it possible to detect them easily (criteria 1). But if multiple outliers are present within \(C(x)\), they might have overly increased the cluster expression interval \(R_{C(x)}\) by inclusion of too many not fitting genes and are not evident in respect to criteria 1. However, because we know that the abundance of gene expression strengths in small RNA-Seq datasets strictly follow a negative exponential distribution (see Section 4.4.1), we can identify expression clusters that include more genes than would be expected from their mean (criteria 2). Outliers are then removed by the double median absolute deviation (MAD) approach, which is is an extension of the normal MAD approach for skewed distributions . After outliers are removed the gene expression clusters can be corrected and are recalculated for every gene. If a single gene \(x\) has too few count values left in any of the investigated conditions, no meaningful assertion can be made about its mean \(\mu_x\) and it should be removed entirely from the differential expression analysis. However, which exact threshold should be used for an exclusion is highly debatable and should always be chosen depending on the general experimental setup. We would recommend to have at least four valid count values per condition for any gene \(x\).

### Hypothesis testing

After the outlier detection and the gene expression cluster correction we suppose that we have good estimators \(\hat{P}_x\) and \(\hat{B}_x\) for every gene \(x\) for two biological conditions I and II. That means we can now model every gene \(x\) as a gamma distributed random variable, depending on its condition: \(\Gamma(\hat{P}_{x_{I}}, \hat{B}_{x_{I}})\) and \(\Gamma(\hat{P}_{x_{II}}, \hat{B}_{x_{II}})\). To evaluate if any gene \(x\) is differentially expressed, *i.e.*, the alternative hypothesis \(\mu_{x_{I}} \ne \mu_{x_{II}}\) is more likely than the null hypothesis \(\mu_{x_{I}} = \mu_{x_{II}}\), we can calculate the probability that the value \(\mu_{x_{I}}\) was drawn from the distribution spanned by \(\Gamma(\hat{P}_{x_{II}}, \hat{B}_{x_{II}})\) (and vice versa). In others words, we calculate the likelihood of drawing the mean count value of gene \(x\) from one condition given the approximated count distribution for the same gene under the second condition. To obtain the respective p-Value we perform three simple steps (see Figure [fig.merde:ablauf]):

1.  Standardization of the gamma distribution spanned by \(\Gamma(\hat{P}_{x_{I}}, \hat{B}_{x_{I}})\) by normalizing the scale parameter \(\hat{B}_{x_{I}}\) with \(\sigma_{x_{I}}\) to obtain \(\Gamma(\hat{P}_{x_{I}}, \dfrac{\hat{B}_{x_{I}}}{\sigma_{x_{I}}})\).

2.  Shift of the standardized gamma distribution by \(\beta\) to obtain \(\Gamma(\hat{P}_{x_{I}}, \beta, \dfrac{\hat{B}_{x_{I}}}{\sigma_{x_{I}}})\), because negative values are not well defined for gamma distributions (see Section 4.4.2, Model description).

3.  Calculation of the p-Value by:

    \[\text{p-Value} = 1 - \int_{-z_{x_{II}}+\beta}^{z_{x_{II}}+\beta} \Gamma(\hat{P}_{x_{I}}, \beta, \dfrac{\hat{B}_{x_{I}}}{\sigma_{x_{I}}}) dz_{x_{II}}\]

    with the z-Values being

    \[z_{x_{II}} = \dfrac{|\mu_{x_{II}} - \mu_{x_{I}}|}{\sigma_{x_{II}}(\sqrt{N_{x_{II}}})^{-1}}\]

    and \(N_{x_{II}}\) being the number of count values of gene \(x\) under condition II.

![image](./figs/merde/ablauf) [fig.merde:ablauf]

There is one exception to this procedure of hypothesis testing. If any given gene \(x\) has too few count values in its corresponding expression cluster \(C(x)\) to calculate good estimators \(\hat{P}_x\) and \(\hat{B}_x\) for its shape and scale parameter, a different approach is used to calculate the respective p-Value. Due to the negative exponential distribution of the expressions strengths in small RNA-Seq data sets, only few miRNA genes with exceptional strong expression are affected by that. As suggest by Qin *et al.* , a cubic root transformation is applied to the count values of the concerning gene \(x\), but instead of a normal t-Test to calculate the p-Value we use Welsch’s t-Test, because we cannot assume an equal variance in the count values of both biological conditions or equal sample sizes . Again, defining how many count values are sufficient to estimate \(\hat{P}_x\) and \(\hat{B}_x\) is debatable, but we recommend to have at least 20 valid count values for any gene \(x\).
