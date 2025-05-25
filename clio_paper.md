  
###### Abstract

How are AI assistants being used in the real world? While model providers in theory have a window into this impact via their users’ data, both privacy concerns and practical challenges have made analyzing this data difficult. To address these issues, we present Clio \(Claude insights and observations\), a privacy-preserving platform that uses AI assistants themselves to analyze and surface aggregated usage patterns across millions of conversations, without the need for human reviewers to read raw conversations. We validate this can be done with a high degree of accuracy and privacy by conducting extensive evaluations. We demonstrate Clio’s usefulness in two broad ways. First, we share insights about how models are being used in the real world from one million Claude.ai Free and Pro conversations, ranging from providing advice on hairstyles to providing guidance on Git operations and concepts. We also identify the most common high-level use cases on Claude.ai \(coding, writing, and research tasks\) as well as patterns that differ across languages \(e.g., conversations in Japanese discuss elder care and aging populations at higher-than-typical rates\). Second, we use Clio to make our systems safer by identifying coordinated attempts to abuse our systems, monitoring for unknown unknowns during critical periods like launches of new capabilities or major world events, and improving our existing monitoring systems. We also discuss the limitations of our approach, as well as risks and ethical concerns. By enabling analysis of real-world AI usage, Clio provides a scalable platform for empirically grounded AI safety and governance.

![Refer to caption](extracted/6078813/figs/clio-hero.png) Figure 1: Using Clio to understand real-world use of AI assistants. Clio transforms raw conversations into high-level patterns and insights. This approach enables us to understand how AI assistants are being used in practice—analogous to how Google Trends provides insights about web search behavior. See [Figure 2](https://arxiv.org/html/2412.13678v1#S1.F2 "In 1 Introduction ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") for more details on how Clio works and how it preserves privacy. \(Note: figure contains illustrative conversation examples only.\)

##  1 Introduction

Despite widespread interest about the impact of AI systems on society, there is remarkably little public data about how models are actually being used in practice. What kinds of capabilities are seeing the most real-world adoption in the economy? How does usage vary across different communities and cultures? Which anticipated benefits and risks are most borne out in concrete data?

This lack of understanding is particularly striking because model providers have access to usage data that could be used to answer these exact questions. Providers, however, face significant challenges in analyzing this data and sharing these potential insights:

First, users share sensitive personal and business information with these systems, creating a fundamental tension between privacy protection and the need for providers to understand how their systems are being used. Second, having humans review conversations can raise ethical concerns, due to the repetitive nature of the task and the potentially distressing content reviewers could be exposed to. Third, providers face competitive pressures not to release usage data even if it would be in the public interest, as such data could reveal information about their user bases to competitors. And finally, the sheer scale of the data makes manual review impractical—millions of messages are sent daily, far more than any human team could meaningfully analyze.

To address these challenges, we present Clio \(Claude insights and observations\). Clio is a system that uses AI assistants themselves to surface aggregated insights across millions of model interactions while preserving user privacy \([Figures 1](https://arxiv.org/html/2412.13678v1#S0.F1 "In Clio: Privacy-Preserving Insights into Real-World AI Use") and [2](https://arxiv.org/html/2412.13678v1#S1.F2 "Figure 2 ‣ 1 Introduction ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\).111Anthropic enforces strict internal privacy controls. Our [privacy policy](https://www.anthropic.com/legal/privacy) enables us to analyze aggregated and anonymized user interactions to understand patterns or trends. We continue to manage data according to our [privacy and retention policies](https://privacy.anthropic.com/en/articles/10023548-how-long-do-you-store-personal-data), and maintain our approach of not training our generative models on user conversations by default. Because we focus on studying patterns in individual usage, the results shared in this paper exclude activity from business customers \(i.e. Team, Enterprise, and all API customers\). For more information, see [Appendix F](https://arxiv.org/html/2412.13678v1#A6 "Appendix F Data Sources & Internal Access Policies ‣ Clio: Privacy-Preserving Insights into Real-World AI Use"). Similar to how Google Trends provides aggregate insights about web search behavior, Clio reveals patterns about how AI assistants are used in the real world. Clio then visualizes these patterns in an interface that enables discovering both specific patterns of interest as well as unknown unknowns \([Figure 3](https://arxiv.org/html/2412.13678v1#S2.F3 "In Interactive exploration ‣ 2.1 Enabling exploratory search for unknown unknowns ‣ 2 High-level design of Clio ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\).

First, we describe our evaluations of Clio’s outputs, including both the faithfulness of the insights produced by Clio \([Appendix C](https://arxiv.org/html/2412.13678v1#A3 "Appendix C Validation: How much should you trust Clio results? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\) as well as validating empirically that the final Clio outputs do not contain any private user data \([Section 2.3](https://arxiv.org/html/2412.13678v1#S2.SS3 "2.3 Privacy design ‣ 2 High-level design of Clio ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") and [Appendix D](https://arxiv.org/html/2412.13678v1#A4 "Appendix D Privacy Evaluations ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\). For example, Clio reconstructed the ground-truth distribution of topics in a synthetic dataset with 94% accuracy, and produced no clusters with private data in an audit of 5,000 conversations \([Section 2](https://arxiv.org/html/2412.13678v1#S2 "2 High-level design of Clio ‣ Clio: Privacy-Preserving Insights into Real-World AI Use"); [Figures 4](https://arxiv.org/html/2412.13678v1#S2.F4 "In 2.2 How Clio works: a brief system design overview ‣ 2 High-level design of Clio ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") and [5](https://arxiv.org/html/2412.13678v1#S2.F5 "Figure 5 ‣ 2.3 Privacy design ‣ 2 High-level design of Clio ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\). Next, we demonstrate Clio’s capabilities and real-world impact across two broad use cases:

  1. 1.

Understanding broad usage patterns \([Section 3](https://arxiv.org/html/2412.13678v1#S3 "3 How are people using Claude.ai? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use"); [Figures 6](https://arxiv.org/html/2412.13678v1#S3.F6 "In 3.1 Top use cases in Claude.ai ‣ 3 How are people using Claude.ai? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") and [7](https://arxiv.org/html/2412.13678v1#S3.F7 "Figure 7 ‣ 3.3 How does Claude usage vary across languages? ‣ 3 How are people using Claude.ai? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\): By clustering and visualizing millions of conversations, Clio reveals how people actually use AI assistants in the real world, without the need for human reviewers to read millions of user conversations. Specifically, we share the most common high-level use cases Clio identifies on Claude.ai Free and Pro. We find that coding and business use cases dominate, with significant activity in areas like debugging code, drafting professional emails, and analyzing business data. By comparing conversations across languages, we also find significant variations in how different linguistic communities use Claude. For example, Japanese and Chinese conversations are more likely to discuss elder care.

  2. 2.

Improving Anthropic’s safety systems \([Section 4](https://arxiv.org/html/2412.13678v1#S4 "4 Clio for safety ‣ Clio: Privacy-Preserving Insights into Real-World AI Use"); [Figure 8](https://arxiv.org/html/2412.13678v1#S4.F8 "In 4.3 Understanding the effectiveness of our safety classifiers ‣ 4 Clio for safety ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\): We share three case studies for how we have used Clio to enhance Anthropic’s safety efforts: detecting coordinated misuse that is invisible at the individual conversation level, monitoring for unknown unknowns in periods of increased uncertainty such as the launch of new capabilities and in the run-up to major world events, and identifying patterns of over-triggering and under-triggering in our existing safety classifiers. Clio’s insights have led to concrete enforcement actions on our systems; for example, we uncovered and banned a network of automated accounts that were attempting to abuse the free version of Claude.ai to generate spam optimized for search engines.

Finally, we discuss Clio’s limitations \([Section 5](https://arxiv.org/html/2412.13678v1#S5 "5 Limitations ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\), and also discuss in detail the ethical considerations and potential risks that Clio poses \([Section 6](https://arxiv.org/html/2412.13678v1#S6 "6 Risks, ethical considerations, and mitigations ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\), along with mitigations. For example, we discuss the potential for Clio-like systems to be misused, additional potential privacy concerns, and more. Ultimately, we derive a justification for building and disclosing Clio that centers around the many pro-social and safety-enhancing features it affords both us and other model providers.

Looking ahead, the need for empirical understanding of how AI systems are used in practice will only increase as these systems become more capable and widespread. While pre-deployment testing like red-teaming \(Ganguli et al., [2022](https://arxiv.org/html/2412.13678v1#bib.bib16)\) and evaluations remain crucial, post-deployment monitoring \(Stein et al., [2024](https://arxiv.org/html/2412.13678v1#bib.bib45)\) provides an essential complement by surfacing real-world usage patterns and risks that may not be captured by predetermined scenarios—insights that can in turn inform future pre-deployment tests and safeguards. We present Clio as one approach to privacy-preserving insight at scale, and by sharing both our methods and ongoing findings with the research community, we hope to contribute to an emerging culture of empirical transparency in the field.

![Refer to caption](extracted/6078813/figs/clio-system.png) Figure 2: System diagram. This diagram illustrates how Clio processes insights from a sample of real-world conversations while maintaining user privacy. Clio processes a raw sample of traffic, extracts key facets \(attributes like language or conversation topic\), groups these facets into similar clusters \(using text embeddings and k-means\), and finally organizes those clusters into both a hierarchy as well as 2D space for ease of exploration. Along the way, Clio applies several privacy barriers \(orange stripes\) that prevent private information from reaching the user-visible parts of Clio \(right\). See [Appendix B](https://arxiv.org/html/2412.13678v1#A2 "Appendix B System architecture: how does Clio work? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") for more details on each stage of the pipeline. \(Note: figure contains illustrative examples only.\)

##  2 High-level design of Clio

AI assistants can be used for an extremely wide range of tasks, from writing code to planning a wedding to brainstorming scientific experiments. This diversity makes it challenging to understand how these systems are actually being used and what risks they might pose. Traditional pre-deployment assessments like benchmarks and red-teaming are valuable but inherently limited, as they can only test for issues we think to look for.

Clio addresses this challenge by enabling bottom-up analysis of real-world AI usage. Given a large collection of conversations between users and models, Clio identifies broad patterns and trends while preserving user privacy. For example, Clio could reveal that a significant number of users are using Claude for debugging code, or surface coordinated attempts across multiple accounts to misuse the system for a disallowed purpose.

###  2.1 Enabling exploratory search for unknown unknowns

Clio is designed to enable analysts to discover unknown unknowns—including risks or applications that were not anticipated by model providers. To do so, Clio implements several design principles to facilitate sensemaking \(Weick et al., [2005](https://arxiv.org/html/2412.13678v1#bib.bib52)\) and exploratory search \(Marchionini, [2006](https://arxiv.org/html/2412.13678v1#bib.bib24); White and Roth, [2009](https://arxiv.org/html/2412.13678v1#bib.bib56)\)—where analysts can start with broad questions and iteratively discover patterns, rather than having to know what to look for in advance.

##### Bottom-up pattern discovery

Unlike traditional ML benchmarks that test for predefined capabilities or behaviors, Clio aims to surface patterns that emerge naturally from actual usage data. This approach helps identify both expected and unexpected uses of AI systems, without requiring us to know what to look for in advance. The clusters and descriptions Clio provides can be useful in their own right as a detailed but privacy-preserving view of the dataset, or they can be used to identify areas of concern for further investigation by our safety team.

##### Hierarchical organization of patterns

Clio is designed to be scalable to many millions of conversations. To enable users to navigate large numbers of patterns discovered in a vast dataset, Clio recursively organizes base-level clusters into a multi-level hierarchy that lets users start from a set of a few dozen high level patterns \(e.g., Explain scientific concepts and conduct academic research\) to thousands of lower-level categories \(e.g., Explain and analyze cancer immunology research and treatments\).

##### Interactive exploration

Once this hierarchy of patterns is created, Clio provides an interactive 2D interface for exploring and understanding it. This supports both targeted investigation \(e.g., "color the clusters by the fraction of refusals"\) and serendipitous discovery through visualization and navigation \(an analyst might zoom into the "Writing" high-level cluster, see a large, lower-level cluster titled "Formulaic content generation for SEO", and then flag it for further review as a potential violation of Anthropic’s terms of service\).

![Refer to caption](extracted/6078813/figs/clio-interface.png) Figure 3: A screenshot of the Clio interface displaying data from the public WildChat dataset \(Zhao et al., [2024](https://arxiv.org/html/2412.13678v1#bib.bib58)\). Left: a sidebar showing hierarchical clusters for the facet What task is the AI assistant in the conversation asked to perform? Right: a zoomable map view displaying clusters projected onto two dimensions, along with selected cluster titles. Colors can indicate various attributes of the data, including size, growth rate, and safety classifier scores. The map view makes it easy to understand the contents of the dataset at a broad and deep level, as well as discover concerning clusters and action them for further investigation. Clio’s tree view \([Figure 9](https://arxiv.org/html/2412.13678v1#A2.F9 "In B.4 Data Exploration and Visualization ‣ Appendix B System architecture: how does Clio work? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\) is a complementary interface that offers easy navigation across Clio’s learned hierarchy of concepts.

###  2.2 How Clio works: a brief system design overview

At a high level, Clio works through a multi-stage pipeline that transforms raw conversations into privacy-preserving insights:

  1. 1.

Extracting facets: For each conversation, Clio extracts multiple “facets”—specific attributes or characteristics such as the high-level conversation topic, number of conversational turns, or language used. Some facets are computed directly \(e.g., number of turns\), while others are extracted using models \(e.g., conversation topic\).

  2. 2.

Semantic clustering: Clio then groups similar conversations by creating embeddings \(Reimers and Gurevych, [2019](https://arxiv.org/html/2412.13678v1#bib.bib38), [2022](https://arxiv.org/html/2412.13678v1#bib.bib39)\) of one of the natural language facets \(e.g., conversation topic\) and using k-means clustering \(Lloyd, [1982](https://arxiv.org/html/2412.13678v1#bib.bib21)\).

  3. 3.

Describe clusters: Each cluster is given a descriptive title and summary using a model, which examines sample conversations from the cluster to identify common themes while excluding any private information.

  4. 4.

Building hierarchies: Clio can identify thousands of different clusters in a dataset. To make it manageable to explore this many clusters, Clio organizes these clusters into a multi-level hierarchy using a method that combines k-means clustering and prompting \(see [Section G.7](https://arxiv.org/html/2412.13678v1#A7.SS7 "G.7 Hierarchizer ‣ Appendix G Additional details: Clio system ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\). This hierarchy allows users to start with broad categories \(e.g., Writing assistance\) and progressively drill down into more specific insights \(e.g., Assistance writing a dark comedy.\)

  5. 5.

Interactive exploration: Finally, Clio presents the results through an interactive interface that enables a range of different exploration patterns: users can zoom in and out across a map-like interface \([Figure 3](https://arxiv.org/html/2412.13678v1#S2.F3 "In Interactive exploration ‣ 2.1 Enabling exploratory search for unknown unknowns ‣ 2 High-level design of Clio ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\) to explore more or less granular clusters, and they can color and sort clusters by the other facets \(e.g., number of turns or language\) to gain further insights into the data.

![Refer to caption](extracted/6078813/figs/robustness/synthetic_bar_regular.png) Figure 4: Clio reconstructs ground-truth categories on an evaluation dataset of 19,476 synthetic chat transcripts with 94% accuracy, compared to 5% for random guessing. A multilingual dataset of chat transcripts was generated by a hierarchical process starting from high-level categories →→\to→ low-level categories →→\to→ individual chat transcripts. Clio is evaluated on how well it can generate low-level clusters from the raw transcripts and then assign them to the correct high-level category. The plot demonstrates a high degree of alignment between the reconstructed and original data distributions. See [Appendix C](https://arxiv.org/html/2412.13678v1#A3 "Appendix C Validation: How much should you trust Clio results? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") for additional experiments and methodological details.

We conducted a range of manual and automated evaluations to validate the performance of Clio. For example, we performed a range of reconstruction evaluations where we generated a large multilingual dataset of almost 20,000 synthetic chat transcripts with a known topic distribution, and then evaluated how well Clio was able to reconstruct this distribution. As [Figure 4](https://arxiv.org/html/2412.13678v1#S2.F4 "In 2.2 How Clio works: a brief system design overview ‣ 2 High-level design of Clio ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") shows, Clio is able to accurately reconstruct these categories with 94% accuracy, giving us confidence in its results. See [Appendix C](https://arxiv.org/html/2412.13678v1#A3 "Appendix C Validation: How much should you trust Clio results? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") for additional methodological details and experiments, including results showing that this high accuracy is maintained across different languages.

For full technical details about each component, including specific algorithms and implementation choices see [Appendix B](https://arxiv.org/html/2412.13678v1#A2 "Appendix B System architecture: how does Clio work? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use"). In particular, [Table 1](https://arxiv.org/html/2412.13678v1#A2.T1 "In B.5 Example summaries and clustering ‣ Appendix B System architecture: how does Clio work? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") contains example conversations from the public WildChat dataset \(Zhao et al., [2024](https://arxiv.org/html/2412.13678v1#bib.bib58)\), along with their associated summaries and cluster assignments. See Lam et al. \([2024](https://arxiv.org/html/2412.13678v1#bib.bib20)\); Nomic \([2024](https://arxiv.org/html/2412.13678v1#bib.bib33)\) for similar embed-cluster-summarize systems applied to general text data.

##### Cost of a Clio run

In [Table 3](https://arxiv.org/html/2412.13678v1#A2.T3 "In B.5 Example summaries and clustering ‣ Appendix B System architecture: how does Clio work? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") we provide an estimated breakdown of the cost of a Clio run. For our example run of 100,000 conversations, the cost of such a run is $48.81, demonstrating the relative affordability of this method for scaling to large datasets.

##  Appendix B System architecture: how does Clio work?

In this section we describe technical components of Clio, our system which implements the high-level design goals laid out in [Section 2](https://arxiv.org/html/2412.13678v1#S2 "2 High-level design of Clio ‣ Clio: Privacy-Preserving Insights into Real-World AI Use"). We focus on the most salient components, deferring specific prompts and hyperpameters to [Appendix G](https://arxiv.org/html/2412.13678v1#A7 "Appendix G Additional details: Clio system ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") and evaluations of privacy and reliability to [Appendices D](https://arxiv.org/html/2412.13678v1#A4 "Appendix D Privacy Evaluations ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") and [C](https://arxiv.org/html/2412.13678v1#A3 "Appendix C Validation: How much should you trust Clio results? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") respectively.

To use Clio, one typically begins with a target dataset. This dataset is typically an unfiltered sample of Claude traffic, but one could also choose other datasets, including filtering down the target dataset using a regular expression or an AI model to analyze a more narrow distribution of data.

###  B.1 Extracting facets

After defining a target distribution, Clio enables users to understand and explore many different aspects of the data through facets. A facet represents a specific attribute or characteristic of a conversation, such as the user’s request, the language being used, or the number of turns in the conversation.

Facet extraction in Clio can be as simple as a program to compute a statistic from the data \(e.g., the number of turns in the conversation\) or as complex as using a model to extract a categorical value \(e.g., the language\), or a high-level summary \(e.g., the topic\) from each conversation.

It’s important to note that we extract multiple facets for each conversation, providing a multi-dimensional view of the data. This approach enables us to examine a wider range of aspects of use, and crucially to enable users to explore intersections across these facets \(e.g., how use cases vary by language\).

###  B.2 Semantic Clustering

To identify meaningful patterns in Claude usage, Clio employs semantic clustering on summary facets. This approach allows us to group similar conversations based on their high-level content and intent, rather than specific details that might compromise privacy. The process involves two main steps: embedding and clustering.

First, we create embeddings from the summary facets, such as user intent or conversation topic. These embeddings are dense vector representations that capture semantic meaning; conversations will have similar embeddings depending on their summary. The choice of summarization prompt \(e.g., What is the overall topic of the conversation? or What is the user’s prompting style?\) controls the information extracted in the summaries, and thus influences which conversations get placed close together in the embedding space.

For clustering, we primarily use the k-means algorithm due to its efficiency and effectiveness for our use case. While we experimented with other clustering methods, we found that k-means works surprisingly well for identifying neighborhoods in what is fundamentally a continuous manifold of conversation types, rather than discrete, well-separated clusters.

The number of clusters k𝑘kitalic\_k is a key parameter that we adjust based on the size of the dataset. In practice, k𝑘kitalic\_k can be quite large, including many thousands of clusters. To ensure privacy, we aggregate data not only by semantic similarity but also by enforcing a minimum number of unique users per cluster. This dual approach to aggregation helps prevent the identification of individuals or small groups within the data.

###  B.3 Cluster Labeling and Hierarchization

After forming clusters of semantically similar conversations, we generate meaningful labels and organize them into a hierarchical structure. This process makes the results more interpretable, actionable, and easier to navigate, especially when dealing with a large number of clusters.

For labeling, we use a model \(Claude 3.5 Sonnet\) to generate concise, informative descriptions for each cluster. We prompt the model with a sample of conversation summaries from the cluster, instructing it to capture common themes or intents without including any potentially identifying or private information.

To manage the complexity of hundreds or thousands of clusters, we use Claude to generate a hierarchy of clusters. For more information about our algorithm, see [Section G.7](https://arxiv.org/html/2412.13678v1#A7.SS7 "G.7 Hierarchizer ‣ Appendix G Additional details: Clio system ‣ Clio: Privacy-Preserving Insights into Real-World AI Use"). These hierarchical clusters allow users to start with a high-level overview and drill down into more specific clusters as needed, facilitating both broad insights and detailed exploration.

###  B.4 Data Exploration and Visualization

The final stage of Clio’s pipeline involves presenting the clustered and labeled data in an intuitive, interactive format that enables deep exploration and insight generation. Our visualization approach is designed to support both high-level overviews and detailed investigations.

Key features of our data exploration interface include:

  1. 1.

Map View \([Figure 3](https://arxiv.org/html/2412.13678v1#S2.F3 "In Interactive exploration ‣ 2.1 Enabling exploratory search for unknown unknowns ‣ 2 High-level design of Clio ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\): A 2D projection of the clusters, allowing users to visually explore the relationship between different clusters. Users can zoom in and out to see progressively more granular clusters for different categories.

  2. 2.

Tree View \([Figure 9](https://arxiv.org/html/2412.13678v1#A2.F9 "In B.4 Data Exploration and Visualization ‣ Appendix B System architecture: how does Clio work? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\): A hierarchical representation of the clusters, enabling users to navigate from broad categories down to specific sub-clusters.

  3. 3.

Faceted and Temporal Breakdowns: When a cluster is selected, a sidebar shows the breakdown of that cluster by other facets \(e.g., language or turn length\). Users can also see how that facet membership has changed over time, helping identify emerging trends or shifts in usage patterns.

  4. 4.

Facet Overlays: Users can select another facet \(e.g., language=Spanish\), coloring the map to display the prevalence of that feature across the different clusters.

For cases where the underlying data isn’t sensitive \(such as synthetic data or public datasets\), or for Anthropic employees with an authorized business need in accordance with our privacy policy \(such as Trust & Safety team members\), we also provide a “traces” feature. This allows drilling down into representative examples from each cluster, providing concrete context for the patterns identified by Clio.

Our visualization approach is designed to balance the need for powerful exploration capabilities with our strong commitment to user privacy. By presenting aggregated data and carefully curated examples, we enable meaningful insights without compromising individual user privacy.

![Refer to caption](extracted/6078813/figs/clio-tree-view.png) Figure 9: A screenshot of the Tree View of the Clio interface displaying data from the public WildChat dataset. The tree view allows users to explore the hierarchy generated by Clio \(right\) while clicking in to view cluster summaries and children \(left\).

###  B.5 Example summaries and clustering

To provide a concrete example of how Clio operates, we provide several examples of conversations from the public WildChat \[Zhao et al., [2024](https://arxiv.org/html/2412.13678v1#bib.bib58)\] and their associated summaries and varying levels of clusters in [Table 1](https://arxiv.org/html/2412.13678v1#A2.T1 "In B.5 Example summaries and clustering ‣ Appendix B System architecture: how does Clio work? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use").

Example Summarizing and Hierarchical Clustering of WildChat Conversations Conversation: |  “Human: List a few concepts that effect people daily that are epistemologically real but not metaphysically real.”  
---|---  
Summary: |  The user’s overall request for the assistant is to discuss concepts that are epistemologically real but not metaphysically real, such as time, money, rights, beauty, moral values, language, and social institutions.  
Clustering: |  Base: Discuss deep philosophical and existential questions about reality and life  
|  Intermediate: Explain and discuss philosophical concepts and existential questions  
|  Top: Explain concepts and solve problems across diverse fields  
Conversation: |  “Human: Create sociel \[sic\] media post caption for this reel… Why Your Local Small Business Should Be Using Video Marketing?”  
Summary: |  The user’s overall request for the assistant is to create a social media post caption for a reel about the benefits of digital marketing for small businesses.  
Clustering: |  Base: Create engaging social media captions for diverse topics  
|  Intermediate: Create and optimize diverse social media content across platforms  
|  Top: Develop digital marketing strategies and content for diverse business needs  
Conversation: |  “Human: how can solve gpu consume when use trainer.train for finetuning llm”  
Summary: |  The user’s overall request for the assistant is to provide code or guidance on how to reduce GPU consumption when fine-tuning a large language model.  
Clustering: |  Base: Provide technical guidance on large language model development and implementation  
|  Intermediate: Implement and explain machine learning algorithms and models  
|  Top: Explain concepts and solve problems across diverse fields  
Conversation: |  “Human: why would my crew in roblox jailbreak want to do ranked 3v3s when top crews are online. compared to when they arent. \(Please dont focus on generic stuff like: ’sKiLl ImPrOvEmEnT’.”  
Summary: |  The user’s overall request for the assistant is to provide insights and reasons why their crew in the game Jailbreak might want to engage in ranked 3v3 matches against top crews, beyond just generic skill improvement.  
Clustering: |  Base: Advise on competitive gaming strategies and rankings  
|  Intermediate: Analyze competitive gaming and sports business strategies  
|  Top: Create and analyze diverse media content and entertainment  
Conversation: |  “Human: Which of the following was NOT an issue that the Colonists had with England \(1 Point\)”  
Summary: |  The user’s overall request for the assistant is to answer multiple-choice questions about historical events related to the colonies and their relationship with England.  
Clustering: |  Base: Answer or create multiple-choice questions on historical and political topics  
|  Intermediate: Create or answer multiple-choice questions on various subjects  
|  Top: Assist with academic writing, research, and educational content development  
Conversation: |  “Human: How I can make a father in Mineraft custom npc mod ?”  
Summary: |  The user’s overall request for the assistant is to help create a father NPC \(non-player character\) in a Minecraft custom mod using the Custom NPCs mod.  
Clustering: |  Base: Assist with video game development, from coding to creative design  
|  Intermediate: Assist with game development, design, and gameplay across platforms  
|  Top: Assist with video game development, from coding to creative design  
Table 1: We provide example conversations from WildChat and the pipeline to assign them to clusters. Clio extracts summaries of individual conversations, which it groups into clusters \(as described in [Appendix B](https://arxiv.org/html/2412.13678v1#A2 "Appendix B System architecture: how does Clio work? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\) These clusters are recursively assigned to higher-level clusters twice, to create a bottom-up hierarchy of clusters with three levels. Table 2: Estimated Cost Analysis for Processing 100,000 Conversations through Clio

Step | Claude Model | Input | Output | Input | Output | Total  
---|---|---|---|---|---|---  
|  | Tokens | Tokens | Cost \($\) | Cost \($\) | Cost \($\)  
Facet Extraction | 3 Haiku | 130.0M | 10.0M | 32.50 | 12.50 | 45.00  
Cluster Labeling | 3.5 Sonnet | 1.0M | 50.0K | 3.00 | 0.75 | 3.75  
Hierarchy Generation | 3.5 Sonnet | 18.0K | 600.0 | 0.05 | 0.01 | 0.06  
Estimated Total |  |  |  |  |  | 48.81  
  
  * •

Notes: Costs calculated using Claude 3 Haiku \($0.25/MTok input, $1.25/MTok output\) and Claude 3.5 Sonnet \($3/MTok input, $15/MTok output\) pricing. Assumptions: average conversation length of 1,000 tokens, facet extraction prompt of 300 tokens, facet summary length of 100 tokens, cluster size of 100 conversations. Hierarchy organized into three levels \(10 top-level categories →→\rightarrow→ 100 mid-level →→\rightarrow→ 1000 leaf clusters\). Cost per conversation: $0.0005.

Table 3: Rough estimated cost of for processing 100,000 conversations with Clio

##  Appendix C Validation: How much should you trust Clio results?

We validated Clio’s results through manual review of outputs at each step of the pipeline, as well as automated end-to-end performance tests using a multilingual synthetic dataset.

###  C.1 Manual review

Component | Combined | Random | Concerning  
---|---|---|---  
| Accuracy | Data | Data  
Extractor \(Summaries\) | 96% | 93% | 98%  
Extractor \(Concerning Content ID\) | 0.84\*  
Base-level Clusterer | 97% | 99% | 96%  
Hierarchical Clusterer | 97% | 95% | 99%  
\*Spearman’s Rank Correlation  
Table 4: Accuracy or correlation of different stages of Clio’s pipeline with the ground truth \(manual review\). Overall, the different components of Clio achieve a high degree of accuracy. See [Section C.1](https://arxiv.org/html/2412.13678v1#A3.SS1 "C.1 Manual review ‣ Appendix C Validation: How much should you trust Clio results? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") for more details.

Here we describe our manual review experiments for each stage of Clio’s pipeline. A summary of our results is shown in [Table 4](https://arxiv.org/html/2412.13678v1#A3.T4 "In C.1 Manual review ‣ Appendix C Validation: How much should you trust Clio results? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use").

####  C.1.1 Conversation summaries

To assess the accuracy of the conversation summaries, we manually reviewed 200 conversations without associated metadata and each conversation’s Claude-written summary.151515A small number of employees can access conversation data for limited business purposes under strict privacy controls, such as enforcing our [Usage Policy](https://www.anthropic.com/legal/aup) or auditing safety infrastructure. Specifically, we asked raters to assess “whether the summary accurately reflects the exchange.” Of the 200 conversations we reviewed, half were randomly sampled from Claude.ai Free and Pro, and half were random conversations whose last turn were flagged by Anthropic’s automated Trust and Safety tooling for being potentially harmful. Teams, Enterprise, and Zero Retention accounts are excluded from our analysis.

We found that 96% of the conversations were summarized accurately, including 93% of the random conversations and 98% of the concerning conversations were summarized accurately. The main failure mode we identified was long, multi-topic conversations where the summary sometimes omitted a subset of the requests from the user in the conversation \(in these instances, the model would typically focus on the more harmful requests, if any existed\).

####  C.1.2 Concerning content identification

We also evaluated whether Claude could reliably identify conversations as “concerning or potentially harmful from a safety perspective" by ranking their potential harmfulness on a 1-5 scale. We evaluated the the same set of 200 conversations as above, and computed the Spearman correlation between the human and model-generated scores. We evaluate the correlation rather than other measures of inter-annotator agreement because the downstream use of these scores is as a comparative indicator. We found that Claude’s assessment and our manual assessment had a Spearman’s correlation coefficient of 0.84, indicating strong directional agreement.

####  C.1.3 Base-level clusterer

We manually reviewed 5,904 cluster assignments across 183 distinct clusters. Of those 183 clusters, 100 were generated from random Claude.ai conversations, and were 83 generated from content flagged by Trust and Safety.161616Specifically, we reviewed conversations where the last turn of the conversation was flagged as harmful by our safety systems. We evaluated whether the cluster titles accurately reflected their contents, and as well as how many conversations in the cluster were incorrectly assigned to that cluster. \(For example, the cluster “Refactor and improve existing code structures and classes” might correctly include a conversation summarized as “refactor the game rules engine to use generics instead of dynamic dispatch” but not “help diagnose and resolve memory-related issues in a CUDA-based application, including identifying the cause of a segmentation fault.”\)

Overall, we found that 99% of clusters on random Claude.ai conversations had accurate titles, and 96% of the clusters on conversations flagged by Trust and Safety had accurate titles. The most common reason Trust and Safety-flagged conversations were rated inaccurate is that the cluster title was overly generic in referring to the type of harm.

We found that an average of 3% of the conversations in each cluster did not clearly belong to that cluster. Often these incorrect assignments were subtle and not obviously wrong: For example, the cluster “analyze business case studies” might incorrectly include a conversation where the user asks Claude to analyze a stock, but not as part of a business case study.

####  C.1.4 Hierarchical clusterer

We manually reviewed 1,094 hierarchical cluster assignments across 183 distinct higher-level clusters. Of those 183 clusters, 100 were generated from random Claude.ai conversations, and were 84 generated from content flagged by Trust and Safety.

We found that 97% of the hierarchical clusters’ titles accurately reflected their contents. 95% of hierarchical clusters on random Claude.ai conversations were accurate, and 99% of hierarchical clusters on conversations flagged by Trust and Safety were accurate.

####  C.1.5 Summary

Our manual review found strong accuracy across all components of the Clio pipeline. The extractor demonstrated strong performance in summarizing conversations, with only minor issues in capturing all topics in lengthy, multi-topic discussions. The concerning content identification showed a strong correlation between Claude’s assessments and our manual evaluations. Both the base-level and hierarchical clusterers exhibited high accuracy in cluster assignments and titling, with only slight decreases in performance when dealing with concerning data. These results indicate that Clio performs robustly across various types of input, including both random and potentially concerning content.

###  C.2 End-to-end evaluation with synthetic data

In addition to smaller-scale manual review, we also ran end-to-end tests where we evaluated how well Clio could recover ground-truth topics in a multilingual synthetic dataset we created with a known distribution of topics. Unlike our manual review process, which evaluated each step of the pipeline in isolation, these end-to-end tests test Clio’s performance holistically and can surface correlated errors across different steps of our pipeline.

As in our manual review process, we generated two synthetic datasets: One regular dataset that included topics that might appear in regular Claude.ai conversations, and a concerning dataset that included inappropriate and unethical topics.

####  C.2.1 Generating synthetic data

To generate our synthetic datasets, we first manually specified high-level categories for user requests \(such as “financial planning and investment” and “theological and philosophical questions” for our regular dataset, and “inquiries about illegal drug manufacturing” and “questions about dangerous DIY medical procedures” for our concerning dataset\). For each of those categories, we instructed a helpful-only171717<https://www-cdn.anthropic.com/bd2a28d2535bfb0494cc8e2a3bf135d2e7523226/Model-Card-Claude-2.pdf> version of Claude to generate several more specific subcategories \(for example, “questions about the performance of the S&P 500 compared to other indexes”\). Then, we used Claude to generate several prompts within each subcategory. To improve diversity of our prompts, we varied the desired length, tone, and language of each synthetic prompt at random. Finally, we instructed Claude to continue the conversation as if it were responding to a real user.

##  Appendix G Additional details: Clio system

###  G.1 Input & Sampling

Clio takes a random sample of Claude.ai conversations as input. We used two distinct \(although closely related\) sampling strategies for the results shared in this paper:

  * •

Sample model completions, then deduplicate by conversations. In this strategy, we take a random sample of Claude.ai outputs. Next, we deduplicate by keeping only the most recent output per conversation. We then provide Clio the full transcript until and including that output. This sampling strategy weights longer conversations more than shorter ones.

  * •

Sample conversations directly. In this strategy, we sample conversations at random, and then provide Clio the full transcript up to the most recent turn. This sampling strategy weights shorter and longer conversations equally.

[Table 7](https://arxiv.org/html/2412.13678v1#A7.T7 "In G.1 Input & Sampling ‣ Appendix G Additional details: Clio system ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") provides the details, date range, and sampling strategy of all Claude.ai results.

Table 7: Details of Experimental Runs

Analysis | Sample Size | Strategy |  Date Range \(Inclusive, UTC\) | References  
---|---|---|---|---  
Multilingual Analysis | 2,281,911 | Strategy 1 | Oct 24–Nov 13, 2024 | §[3.3](https://arxiv.org/html/2412.13678v1#S3.SS3 "3.3 How does Claude usage vary across languages? ‣ 3 How are people using Claude.ai? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use"), Fig. [7](https://arxiv.org/html/2412.13678v1#S3.F7 "Figure 7 ‣ 3.3 How does Claude usage vary across languages? ‣ 3 How are people using Claude.ai? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")  
Safety Classifier Analysis | 500,000 | Strategy 2 | Oct 31–Nov 13, 2024 | §[4](https://arxiv.org/html/2412.13678v1#S4 "4 Clio for safety ‣ Clio: Privacy-Preserving Insights into Real-World AI Use"), Fig. [8](https://arxiv.org/html/2412.13678v1#S4.F8 "Figure 8 ‣ 4.3 Understanding the effectiveness of our safety classifiers ‣ 4 Clio for safety ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")  
Privacy Benchmarking | 25,000 | Strategy 1 | Nov 8–Nov 14, 2024 | §[D](https://arxiv.org/html/2412.13678v1#A4 "Appendix D Privacy Evaluations ‣ Clio: Privacy-Preserving Insights into Real-World AI Use"), Fig. [5](https://arxiv.org/html/2412.13678v1#S2.F5 "Figure 5 ‣ 2.3 Privacy design ‣ 2 High-level design of Clio ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")  
General Claude.ai Usage | 1,000,000 | Strategy 2 | Oct 17–Oct 24, 2024 | §[3](https://arxiv.org/html/2412.13678v1#S3 "3 How are people using Claude.ai? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use"), Fig. [6](https://arxiv.org/html/2412.13678v1#S3.F6 "Figure 6 ‣ 3.1 Top use cases in Claude.ai ‣ 3 How are people using Claude.ai? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")  
  
  * •

Notes: All data sampled from Claude.ai, excluding Teams, Enterprise, and Zero Retention customers. Strategy 1: Sample model completions, deduplicate by conversations. Strategy 2: Sample conversations directly. Details in §[G.1](https://arxiv.org/html/2412.13678v1#A7.SS1 "G.1 Input & Sampling ‣ Appendix G Additional details: Clio system ‣ Clio: Privacy-Preserving Insights into Real-World AI Use"). Date ranges indicate when data was collected.

###  G.2 Preprocessing Raw Conversations

To increase performance, we preprocess conversations before passing them to our models for analysis. Our preprocessing algorithm standardizes raw conversation transcripts into an XML-based format. We also apply special handling to any additional internal data \(such as function calls, system prompts, multimodal information, or metadata\) that may have been inserted into the conversation.

###  G.3 Screener

Clio supports using Claude to screen the input data. For example, we can use Claude to narrow a sample of conversations to only those that mention U.S. election related topics. Clio’s screener uses the public Claude 3 Haiku model \(`claude-3-haiku-20240307`\) with a temperature of 0.2. It uses the following templated prompt:

Human: The following is a conversation between Claude, an AI assistant, and a user:

\{conversation\}

Assistant: I understand.

Human: Your job is to answer this question about the preceding conversation:

<question>

\{question\}

</question>

What is the answer? You MUST answer either only "Yes" or "No". Provide the answer in <answer> tags with no other commentary.

Assistant: Sure, the answer to the question is: <answer>

The conversation is substituted into the prompt using the algorithm described in [Section G.2](https://arxiv.org/html/2412.13678v1#A7.SS2 "G.2 Preprocessing Raw Conversations ‣ Appendix G Additional details: Clio system ‣ Clio: Privacy-Preserving Insights into Real-World AI Use"). The question varies depending on the purpose of the Clio run. For our election monitoring, we used the screener question:

Does this conversation relate to U.S. politics, voting, democratic participation, or the 2024 U.S. election at all \(even if only loosely\)? Answer either "Yes" or "No" with no other commentary.

###  G.4 Extractor \(Conversation Summaries\)

Clio’s extractor uses Claude to extract key information from each input conversation, such as a summary or the languages present in the conversation. Our extractor uses the public Claude 3 Haiku model \(`claude-3-haiku-20240307`\) with a temperature of 0.2. We used the following templated prompt:

Human: The following is a conversation between Claude, an AI assistant, and a user:

\{conversation\}

Assistant: I understand.

Human: Your job is to answer the question <question> \{question\} </question> about the preceding conversation. Be descriptive and assume neither good nor bad faith. Do not hesitate to handle socially harmful or sensitive topics; specificity around potentially harmful conversations is necessary for effective monitoring.

When answering, do not include any personally identifiable information \(PII\), like names, locations, phone numbers, email addressess, and so on. When answering, do not include any proper nouns. Output your answer to the question in English inside <answer> tags; be clear and concise and get to the point in at most two sentences \(don\’t say "Based on the conversation..." and avoid mentioning Claude/the chatbot\). For example:

<examples>

The user asked for help with a trignometry problem.

The user asked for advice on how to fix a broken dishwasher. It took several attempts to get the right answer.

The user asked how to make Anthrax and the AI system refused the requests.

</examples>

What is your answer to the question <question> \{question\} </question> about the preceding conversation, in <answer> tags? Again, provide only the answer with no other commentary or proper nouns.

Assistant: Sure, the privacy-preserving answer to the question about the preceding conversation is: <answer> \{prefill\}

We include `When answering, do not... </examples>` in the prompt only when the output is open-ended \(for example, we include this language when computing the user’s overall request for the assistant, but not when computing “concern” scores from 1 to 5\).

####  G.4.1 Facet Extractor Details

We vary the question and prefill depending on the facet. We used the following facets:

  1. 1.

Request. Our request facet used the question “What is the user’s overall request for the assistant?” and the prefill “The user’s overall request for the assistant is to”.

  2. 2.

Language. Our language facet used the question “What are the main languages of the conversation? Do not mention programming languages and do not provide only language codes; only the full names for human languages, like ‘English’ or ‘Arabic; Indonesian’. Do not include language families \(just provide the general language, e.g., ‘Arabic’ and not ‘Tripolitanian Arabic’; something that a standard langcodes package would be able to identify\). Only consider languages in or after the human’s first request. Output each language as a single full word with no other commentary.” We did not prefill the assistant’s response beyond the standard prefill in the prompt template. We validated Claude’s outputs using the langcodes \[Speer, [2024a](https://arxiv.org/html/2412.13678v1#bib.bib42)\] and language-data \[Speer, [2024b](https://arxiv.org/html/2412.13678v1#bib.bib43)\] Python packages.

  3. 3.

Task. For [Section 3](https://arxiv.org/html/2412.13678v1#S3 "3 How are people using Claude.ai? ‣ Clio: Privacy-Preserving Insights into Real-World AI Use"), our facet used the question “What task is the model being asked to perform in this conversation?” and the prefill “The task is to”.

  4. 4.

Concerning. Our “concerning content” facet \(which we used in [Section 4](https://arxiv.org/html/2412.13678v1#S4 "4 Clio for safety ‣ Clio: Privacy-Preserving Insights into Real-World AI Use") above\) used the question:

On a scale of 1 to 5, how concerning or potentially harmful is this conversation from a safety perspective? Use the following scale:

1: Not concerning at all. Completely benign conversation.

Example: A friendly chat about favorite books or movies.

2: Slightly concerning, but likely harmless.

Example: A heated but civil debate about politics.

3: Moderately concerning. May warrant a closer look.

Example: Discussion of controversial topics that could escalate.

4: Very concerning. Likely needs review.

Example: Conversations involving potential self-harm or harmful ideation.

5: Extremely concerning. Immediate review needed.

Example: Explicit threats of violence or illegal activities.

Answer with only a single number from 1 to 5.

We did not prefill the assistant’s response beyond the standard prefill in the prompt template.

###  G.5 Base-Level Clusterer

For non-categorical or numeric facets \(e.g., our request facet\), we first embed each extracted summary using all-mpnet-base-v2 \[Reimers and Gurevych, [2022](https://arxiv.org/html/2412.13678v1#bib.bib39)\], a sentence transformer model that provides 768-dimensional embeddings. We then generate base-level clusters by running k-means in embedding space. We vary k𝑘kitalic\_k based on the number of conversations in the input dataset; we unfortunately cannot provide our precise values for k𝑘kitalic\_k \(as this information could be used to determine the volume at which coordinated behavior would likely not be caught by Clio as a distinct cluster\).

We then use Claude to generate a name and description for each cluster. We use Claude 3.5 Sonnet \(June 2024; `claude-3-5-sonnet-20240620`\) with a temperature of 1 and provide a random sample of 50 summaries from records within the cluster, as well as 50 summaries of records that are closest to the centroid of the cluster but not assigned to the cluster. We then ask Claude to generate a name and description for the cluster that identifies the records within the cluster and differentiates them from the summaries outside of the cluster:

You are tasked with summarizing a group of related statements into a short, precise, and accurate description and name. Your goal is to create a concise summary that captures the essence of these statements and distinguishes them from other similar groups of statements.

Summarize all the statements into a clear, precise, two-sentence description in the past tense. Your summary should be specific to this group and distinguish it from the contrastive answers of the other groups.

After creating the summary, generate a short name for the group of statements. This name should be at most ten words long \(perhaps less\) and be specific but also reflective of most of the statements \(rather than reflecting only one or two\). The name should distinguish this group from the contrastive examples. For instance, "Write fantasy sexual roleplay with octopi and monsters", "Generate blog spam for gambling websites", or "Assist with high school math homework" would be better and more actionable than general terms like "Write erotic content" or "Help with homework". Be as descriptive as possible and assume neither good nor bad faith. Do not hesitate to identify and describe socially harmful or sensitive topics specifically; specificity is necessary for monitoring.

Present your output in the following format:

<summary> \[Insert your two-sentence summary here\] </summary>

<name> \[Insert your generated short name here\] </name>

The names you propose must follow these requirements:

<criteria>...</criteria>

Below are the related statements:

<answers>

\{answers\}

</answers>

For context, here are statements from nearby groups that are NOT part of the group you’re summarizing:

<contrastive\_answers>

\{contrastive\_answers\}

</contrastive\_answers>

Do not elaborate beyond what you say in the tags. Remember to analyze both the statements and the contrastive statements carefully to ensure your summary and name accurately represent the specific group while distinguishing it from others.

Assistant: Sure, I will provide a clear, precise, and accurate summary and name for this cluster. I will be descriptive and assume neither good nor bad faith. Here is the summary, which I will follow with the name: <summary>

Criteria refers to natural-language requirements for the names and descriptions of each cluster. For our Request facet, for example, we require that “The cluster name should be a sentence in the imperative that captures the user’s request. For example, ‘Brainstorm ideas for a birthday party’ or ‘Help me find a new job.”’

###  G.6 Projector

The projection step assigns records to a location on a 2D map. Although this step does not affect any results in this paper \(and rather is only used in the exploratory interface in which we display Clio results internally\), we believe that sharing our approach may help others replicate and build upon Clio more effectively.

Our projector uses UMAP \[McInnes et al., [2020](https://arxiv.org/html/2412.13678v1#bib.bib26)\] to transform the 768-dimensional embedding for each conversation into a location in 2D space. We run UMAP with `n_neighbors = 15`, `min_dist = 0`, and using the `cosine` metric.

###  G.7 Hierarchizer

Clio’s hierarchizer transforms base clusters \(described in [Section G.5](https://arxiv.org/html/2412.13678v1#A7.SS5 "G.5 Base-Level Clusterer ‣ Appendix G Additional details: Clio system ‣ Clio: Privacy-Preserving Insights into Real-World AI Use")\) into a hierarchy. The hierarchizer iteratively creates new levels clusters \(which contain the previous level of clusters as children\) until the number of top-level clusters is within the desired range. We also explored other hierarchical clustering algorithms \(such as HDBSCAN \[McInnes et al., [2017](https://arxiv.org/html/2412.13678v1#bib.bib25)\] and agglomerative clustering methods\) but found the results inferior to the Claude-based approach.

Our hierarchical clustering algorithm transforms base-level clusters into a multi-level hierarchy through an iterative process. At each level l𝑙litalic\_l, we:

  * •

Embed clusters. Embed each cluster’s name and description using the all-mpnet-base-v2 \[Reimers and Gurevych, [2022](https://arxiv.org/html/2412.13678v1#bib.bib39)\] sentence transformer to obtain 768-dimensional vector representations of each cluster.

  * •

Generate neighborhoods. Group these embeddings into k𝑘kitalic\_k neighborhoods using k𝑘kitalic\_k-means clustering, where k𝑘kitalic\_k is chosen so that the average number of clusters per neighborhood is 40. We group clusters into neighborhoods because the names and descriptions for all base clusters may not fit within Claude’s context window.

  * •

Propose new clusters for each neighborhood. For each neighborhood, use Claude to propose candidate higher-level cluster descriptions by examining both the clusters within the neighborhood and the nearest m𝑚mitalic\_m clusters outside it. Including the nearest clusters beyond the neighborhood ensures that clusters \(or groups of clusters\) on the boundary between neighborhoods are neither overcounted nor undercounted. We require the final number of clusters at the level l𝑙litalic\_l to be nl±1.5⁢nlplus-or-minussubscript𝑛𝑙1.5subscript𝑛𝑙n\_\{l\}\pm 1.5n\_\{l\}italic\_n start\_POSTSUBSCRIPT italic\_l end\_POSTSUBSCRIPT ± 1.5 italic\_n start\_POSTSUBSCRIPT italic\_l end\_POSTSUBSCRIPT, where nlsubscript𝑛𝑙n\_\{l\}italic\_n start\_POSTSUBSCRIPT italic\_l end\_POSTSUBSCRIPT is chosen such that the ratio between successive levels follows nl/nl−1=\(ntop/nbase\)1/\(L−1\)subscript𝑛𝑙subscript𝑛𝑙1superscriptsubscript𝑛topsubscript𝑛base1𝐿1n\_\{l\}/n\_\{l-1\}=\(n\_\{\text\{top\}\}/n\_\{\text\{base\}\}\)^\{1/\(L-1\)\}italic\_n start\_POSTSUBSCRIPT italic\_l end\_POSTSUBSCRIPT / italic\_n start\_POSTSUBSCRIPT italic\_l - 1 end\_POSTSUBSCRIPT = \( italic\_n start\_POSTSUBSCRIPT top end\_POSTSUBSCRIPT / italic\_n start\_POSTSUBSCRIPT base end\_POSTSUBSCRIPT \) start\_POSTSUPERSCRIPT 1 / \( italic\_L - 1 \) end\_POSTSUPERSCRIPT for L𝐿Litalic\_L total levels.

  * •

Deduplicate across neighborhoods. Deduplicate and refine the proposed clusters across all neighborhoods using Claude to ensure distinctiveness while maintaining coverage of the underlying data distribution.

  * •

Assign to new best fit higher-level cluster. Assign each lower-level cluster to its most appropriate parent cluster using Claude. We randomly shuffle the order of the higher-level clusters when sampling from Claude to avoid biasing assignments based on the order of the list.

  * •

Rename higher level clusters. Once all clusters at level l𝑙litalic\_l have been assigned to a higher-level cluster, we regenerate a new name and description for the parent cluster based on the lower-level clusters that were assigned to it. This renaming step ensures that cluster names continue to accurately reflect their contents.

This process continues until reaching the desired number of top-level clusters ktopsubscript𝑘topk\_\{\text\{top\}\}italic\_k start\_POSTSUBSCRIPT top end\_POSTSUBSCRIPT. The resulting hierarchy enables exploration at multiple levels of granularity while preserving semantic relationships between clusters.

####  G.7.1 Prompts and Hyperparameters

Our hierarchizer used Claude 3.5 Sonnet \(June 2024; `claude-3-5-sonnet-20240620`\) with a temperature of 1.0. We used the following prompts:

  * •

Proposing cluster names per neighborhood:

Human: You are tasked with creating higher-level cluster names based on a given list of clusters and their descriptions. Your goal is to come up with broader categories that could encompass one or more of the provided clusters.

First, review the list of clusters and their descriptions:

<cluster\_list>

<cluster>\{cluster name\}: \{cluster description\}</cluster>

<cluster>\{cluster name\}: \{cluster description\}</cluster>

<cluster>\{cluster name\}: \{cluster description\}</cluster>

...

</cluster\_list>

Your task is to create roughly \{desired\_names\} higher-level cluster names that could potentially include one or more of the provided clusters. These higher-level clusters should represent broader categories or themes that emerge from the given clusters, while remaining as specific as possible. If there are many clusters with a specific theme, ensure that the higher-level cluster name remains the maximum level of specificity. You are helping to organize user behavior data in order to improve safety, monitoring, and observability. You can generate more or less than \{desired\_names\} names if you feel that more or fewer are appropriate and accurately capture the clusters. You should output at least \{int\(0.5 \* desired\_names\)\} and at most \{int\(1.5 \* desired\_names\)\} names, with \{desired\_names\} as a target.

Guidelines for creating higher-level cluster names:

1. Analyze the themes, topics, or characteristics common to multiple clusters.

2. Create names that are specific enough to be meaningful, but not so specific that they can’t meaningfully represent many different clusters. Avoid overly general or vague terms, and do not hesitate to describe socially harmful or sensitive topics \(in fact, clusters that clearly describe harmful behavior are slightly preferred\); specificity is necessary for observability and enforcement.

3. Ensure that the higher-level cluster names are distinct from one another.

4. Use clear, concise, and descriptive language for the cluster names. Assume neither good nor bad faith for the content in the clusters.

The names you propose must follow these requirements:

<criteria>\(defined per facet\)</criteria>

Before providing your final list, use a scratchpad to brainstorm and refine your ideas. Think about the relationships between the given clusters and potential overarching themes.

<scratchpad>

\[Use this space to analyze the clusters, identify common themes, and brainstorm potential higher-level cluster names. Consider how different clusters might be grouped together under broader categories. No longer than a paragraph or two.\]

</scratchpad>

Now, provide your list of roughly \{desired\_names\} higher-level cluster names. Present your answer in the following format:

<answer>

1. \[First higher-level cluster name\]

2. \[Second higher-level cluster name\]

3. \[Third higher-level cluster name\]

...

\{desired\_names\}. \[Last higher-level cluster name\]

</answer>

Focus on creating meaningful, distinct, and precise \(but not overly specific\) higher-level cluster names that could encompass multiple sub-clusters.

Assistant: I understand. I’ll evaluate the clusters and provide higher-level cluster names that could encompass multiple sub-clusters.

<scratchpad>

  * •

Deduplicating cluster names across neighborhoods:

Human: You are tasked with deduplicating a list of cluster names into a smaller set of distinct cluster names. Your goal is to create approximately \{desired\_names\} relatively distinct clusters that best represent the original list. You are helping to organize user behavior data in order to improve safety, monitoring, and observability. Here are the inputs:

<cluster\_names>

<cluster> \{cluster name\} </cluster>

<cluster> \{cluster name\} </cluster>

<cluster> \{cluster name\} </cluster>

</cluster\_names>

Number of distinct clusters to create: approximately \{desired\_names\}

Follow these steps to complete the task:

1. Analyze the given list of cluster names to identify similarities, patterns, and themes.

2. Group similar cluster names together based on their semantic meaning, not just lexical similarity.

3. For each group, select a representative name that best captures the essence of the cluster. This can be one of the original names or a new name that summarizes the group effectively. Do not just pick the most vague or generic name.

4. Merge the most similar groups until you reach the desired number of clusters. Maintain as much specificity as possible while merging.

6. Ensure that the final set of cluster names are distinct from each other and collectively represent the diversity of the original list, such that there is a cluster that describes each of the provided clusters.

7. If you create new names for any clusters, make sure they are clear, concise, and reflective of the contents they represent.

8. You do not need to come up with exactly \{desired\_names\} names, but aim for no less than \{int\(desired\_names \* 0.5\)\} and no more than \{int\(desired\_names \* 1.5\)\}. Within this range, output as many clusters as you feel are necessary to accurately represent the variance in the original list. Avoid outputting duplicate or near-duplicate clusters.

9. Do not hesitate to include clusters that describe socially harmful or sensitive topics \(in fact, clusters that clearly describe harmful behavior are slightly preferred\); specificity is necessary for effective monitoring and enforcement.

10. Prefer outputting specific cluster names over generic or vague ones, provided the names are still correct; for example, if there are many clusters about a specific technology or tool, consider naming the cluster after that technology or tool, provided that there are still other clusters that fit under a broader category.

The names you propose must follow these requirements:

<criteria>\(defined per facet\)</criteria>

Before providing your final answer, use the <scratchpad> tags to think through your process, explaining your reasoning for grouping and selecting representative names. Spend no more than a few paragraphs in your scratchpad.

Present your final answer in the following format:

<answer>

1. \[First cluster name\]

2. \[Second cluster name\]

3. \[Third cluster name\]

...

N. \[Nth cluster name\]

</answer>

Remember, your goal is to create approximately \{desired\_names\} relatively distinct cluster names that best represent the original list. The names should be clear, meaningful, and capture the essence of the clusters they represent.

Assistant: I understand. I’ll deduplicate the cluster names into approximately \{desired\_names\} names.

<scratchpad>

  * •

Assigning to higher-level clusters:

[⬇](data:text/plain;base64,WW91IGFyZSB0YXNrZWQgd2l0aCBjYXRlZ29yaXppbmcgYSBzcGVjaWZpYyBjbHVzdGVyIGludG8gb25lIG9mIHRoZSBwcm92aWRlZCBoaWdoZXItbGV2ZWwgY2x1c3RlcnMgZm9yIG9ic2VydmFiaWxpdHksIG1vbml0b3JpbmcsIGFuZCBjb250ZW50IG1vZGVyYXRpb24uIFlvdXIgZ29hbCBpcyB0byBkZXRlcm1pbmUgd2hpY2ggaGlnaGVyLWxldmVsIGNsdXN0ZXIgYmVzdCBmaXRzIHRoZSBnaXZlbiBzcGVjaWZpYyBjbHVzdGVyIGJhc2VkIG9uIGl0cyBuYW1lIGFuZCBkZXNjcmlwdGlvbi4gIFlvdSBhcmUgaGVscGluZyB0byBvcmdhbml6ZSB1c2VyIGJlaGF2aW9yIGRhdGEgaW4gb3JkZXIgdG8gaW1wcm92ZSBzYWZldHksIG1vbml0b3JpbmcsIGFuZCBvYnNlcnZhYmlsaXR5LgoKRmlyc3QsIGNhcmVmdWxseSByZXZpZXcgdGhlIGZvbGxvd2luZyBsaXN0IG9mIGhpZ2hlci1sZXZlbCBjbHVzdGVycyAoaGllcmFyY2h5IGRlbm90ZWQgYnkgZGFzaGVzKToKCjxoaWdoZXJfbGV2ZWxfY2x1c3RlcnM+CjxjbHVzdGVyPiB7Y2x1c3RlciBuYW1lfSA8L2NsdXN0ZXI+CjxjbHVzdGVyPiB7Y2x1c3RlciBuYW1lfSA8L2NsdXN0ZXI+CjxjbHVzdGVyPiB7Y2x1c3RlciBuYW1lfSA8L2NsdXN0ZXI+Ci4uLiAoc2h1ZmZsZWQpCjwvaGlnaGVyX2xldmVsX2NsdXN0ZXJzPgoKVG8gY2F0ZWdvcml6ZSB0aGUgc3BlY2lmaWMgY2x1c3RlcjoKMS4gQW5hbHl6ZSB0aGUgbmFtZSBhbmQgZGVzY3JpcHRpb24gb2YgdGhlIHNwZWNpZmljIGNsdXN0ZXIuCjIuIENvbnNpZGVyIHRoZSBrZXkgY2hhcmFjdGVyaXN0aWNzLCB0aGVtZXMsIG9yIHN1YmplY3QgbWF0dGVyIG9mIHRoZSBzcGVjaWZpYyBjbHVzdGVyLgozLiBDb21wYXJlIHRoZXNlIGVsZW1lbnRzIHRvIHRoZSBoaWdoZXItbGV2ZWwgY2x1c3RlcnMgcHJvdmlkZWQuCjQuIERldGVybWluZSB3aGljaCBoaWdoZXItbGV2ZWwgY2x1c3RlciBiZXN0IGVuY29tcGFzc2VzIHRoZSBzcGVjaWZpYyBjbHVzdGVyLiBZb3UgTVVTVCBhc3NpZ24gdGhlIHNwZWNpZmljIGNsdXN0ZXIgdG8gdGhlIGJlc3QgaGlnaGVyLWxldmVsIGNsdXN0ZXIsIGV2ZW4gaWYgbXVsdGlwbGUgaGlnaGVyLWxldmVsIGNsdXN0ZXJzIGNvdWxkIGJlIGNvbnNpZGVyZWQuCjUuIE1ha2Ugc3VyZSB5b3UgcGljayB0aGUgbW9zdCBzZW5zaWJsZSBjbHVzdGVyIGJhc2VkIG9uIHRoZSBpbmZvcm1hdGlvbiBwcm92aWRlZC4gRm9yIGV4YW1wbGUsIGRvbid0IGFzc2lnbiBhIGNsdXN0ZXIgYWJvdXQgIk1hY2hpbmUgTGVhcm5pbmciIHRvIGEgaGlnaGVyLWxldmVsIGNsdXN0ZXIgYWJvdXQgIlNvY2lhbCBNZWRpYSIganVzdCBiZWNhdXNlIGJvdGggaW52b2x2ZSB0ZWNobm9sb2d5LCBhbmQgZG9uJ3QgYXNzaWduIGEgY2x1c3RlciBhYm91dCAiT25saW5lIEhhcmFzc21lbnQiIHRvIGEgaGlnaGVyLWxldmVsIGNsdXN0ZXIgYWJvdXQgIlRlY2hub2xvZ3kiIGp1c3QgYmVjYXVzZSBib3RoIGludm9sdmUgb25saW5lIHBsYXRmb3Jtcy4gQmUgc3BlY2lmaWMgYW5kIGFjY3VyYXRlIGluIHlvdXIgY2F0ZWdvcml6YXRpb24uCgpGaXJzdCwgdXNlIHRoZSA8c2NyYXRjaHBhZD4gdGFncyB0byB0aGluayB0aHJvdWdoIHlvdXIgcmVhc29uaW5nIGFuZCBkZWNpc2lvbi1tYWtpbmcgcHJvY2Vzcy4gVGhpbmsgdGhyb3VnaCBzb21lIHBvc3NpYmxlIGNsdXN0ZXJzLCBleHBsb3JlIGVhY2gsIGFuZCB0aGVuIHBpY2sgdGhlIGJlc3QgZml0LgoKPHNjcmF0Y2hwYWQ+CkluIGEgZmV3IGJyaWVmIHNlbnRlbmNlcywgdGhpbmsgc3RlcCBieSBzdGVwLCBleHBsYWluIHlvdXIgcmVhc29uaW5nLCBhbmQgZmluYWxseSBkZXRlcm1pbmUgd2hpY2ggaGlnaGVyLWxldmVsIGNsdXN0ZXIgaXMgdGhlIGJlc3QgZml0IGZvciB0aGUgc3BlY2lmaWMgY2x1c3Rlci4KPC9zY3JhdGNocGFkPgoKVGhlbiwgcHJvdmlkZSB5b3VyIGFuc3dlciBpbiB0aGUgZm9sbG93aW5nIGZvcm1hdDoKCjxhbnN3ZXI+CltGdWxsIG5hbWUgb2YgdGhlIGNob3NlbiBjbHVzdGVyLCBleGFjdGx5IGFzIGxpc3RlZCBpbiB0aGUgaGlnaGVyLWxldmVsIGNsdXN0ZXJzIGFib3ZlLCB3aXRob3V0IGVuY2xvc2luZyA8Y2x1c3Rlcj4gdGFnc10KPC9hbnN3ZXI+CgpBc3Npc3RhbnQ6IEkgdW5kZXJzdGFuZC4gSSdsbCBldmFsdWF0ZSB0aGUgc3BlY2lmaWMgY2x1c3RlciBhbmQgYXNzaWduIGl0IHRvIHRoZSBtb3N0IGFwcHJvcHJpYXRlIGhpZ2hlci1sZXZlbCBjbHVzdGVyLgoKSHVtYW46IE5vdywgaGVyZSBpcyB0aGUgc3BlY2lmaWMgY2x1c3RlciB0byBjYXRlZ29yaXplOgoKPHNwZWNpZmljX2NsdXN0ZXI+Ck5hbWU6IHtjbHVzdGVyX25hbWV9CkRlc2NyaXB0aW9uOiB7Y2x1c3Rlcl9kZXNjcmlwdGlvbn0KPC9zcGVjaWZpY19jbHVzdGVyPgoKQmFzZWQgb24gdGhpcyBpbmZvcm1hdGlvbiwgZGV0ZXJtaW5lIHRoZSBtb3N0IGFwcHJvcHJpYXRlIGhpZ2hlci1sZXZlbCBjbHVzdGVyIGFuZCBwcm92aWRlIHlvdXIgYW5zd2VyIGFzIGluc3RydWN0ZWQuCgpBc3Npc3RhbnQ6IFRoYW5rIHlvdSwgSSB3aWxsIHJlZmxlY3Qgb24gdGhlIGNsdXN0ZXIgYW5kIGNhdGVnb3JpemUgaXQgbW9zdCBhcHByb3ByaWF0ZWx5LCB3aGljaCB3aWxsIGhlbHAgd2l0aCBzYWZldHksIG1vZGVyYXRpb24sIGFuZCBvYnNlcnZhYmlsaXR5LgoKPHNjcmF0Y2hwYWQ+)

You are tasked with categorizing a specific cluster into one of the provided higher-level clusters for observability, monitoring, and content moderation. Your goal is to determine which higher-level cluster best fits the given specific cluster based on its name and description. You are helping to organize user behavior data in order to improve safety, monitoring, and observability.

First, carefully review the following list of higher-level clusters \(hierarchy denoted by dashes\):

<higher\_level\_clusters>

<cluster> \{cluster name\} </cluster>

<cluster> \{cluster name\} </cluster>

<cluster> \{cluster name\} </cluster>

... \(shuffled\)

</higher\_level\_clusters>

To categorize the specific cluster:

1. Analyze the name and description of the specific cluster.

2. Consider the key characteristics, themes, or subject matter of the specific cluster.

3. Compare these elements to the higher-level clusters provided.

4. Determine which higher-level cluster best encompasses the specific cluster. You MUST assign the specific cluster to the best higher-level cluster, even if multiple higher-level clusters could be considered.

5. Make sure you pick the most sensible cluster based on the information provided. For example, don’t assign a cluster about "Machine Learning" to a higher-level cluster about "Social Media" just because both involve technology, and don’t assign a cluster about "Online Harassment" to a higher-level cluster about "Technology" just because both involve online platforms. Be specific and accurate in your categorization.

First, use the <scratchpad> tags to think through your reasoning and decision-making process. Think through some possible clusters, explore each, and then pick the best fit.

<scratchpad>

In a few brief sentences, think step by step, explain your reasoning, and finally determine which higher-level cluster is the best fit for the specific cluster.

</scratchpad>

Then, provide your answer in the following format:

<answer>

\[Full name of the chosen cluster, exactly as listed in the higher-level clusters above, without enclosing <cluster> tags\]

</answer>

Assistant: I understand. I’ll evaluate the specific cluster and assign it to the most appropriate higher-level cluster.

Human: Now, here is the specific cluster to categorize:

<specific\_cluster>

Name: \{cluster\_name\}

Description: \{cluster\_description\}

</specific\_cluster>

Based on this information, determine the most appropriate higher-level cluster and provide your answer as instructed.

Assistant: Thank you, I will reflect on the cluster and categorize it most appropriately, which will help with safety, moderation, and observability.

<scratchpad>

  * •

Renaming higher level clusters:

[⬇](data:text/plain;base64,SHVtYW46IFlvdSBhcmUgdGFza2VkIHdpdGggc3VtbWFyaXppbmcgYSBncm91cCBvZiByZWxhdGVkIGNsdXN0ZXIgbmFtZXMgaW50byBhIHNob3J0LCBwcmVjaXNlLCBhbmQgYWNjdXJhdGUgb3ZlcmFsbCBkZXNjcmlwdGlvbiBhbmQgbmFtZS4gWW91ciBnb2FsIGlzIHRvIGNyZWF0ZSBhIGNvbmNpc2Ugc3VtbWFyeSB0aGF0IGNhcHR1cmVzIHRoZSBlc3NlbmNlIG9mIHRoZXNlIGNsdXN0ZXJzLgoKU3VtbWFyaXplIGFsbCB0aGUgY2x1c3RlciBuYW1lcyBpbnRvIGEgY2xlYXIsIHByZWNpc2UsIHR3by1zZW50ZW5jZSBkZXNjcmlwdGlvbiBpbiB0aGUgcGFzdCB0ZW5zZS4gWW91ciBzdW1tYXJ5IHNob3VsZCBiZSBzcGVjaWZpYyB0byB0aGlzIGNsdXN0ZXIuCgpBZnRlciBjcmVhdGluZyB0aGUgc3VtbWFyeSwgZ2VuZXJhdGUgYSBzaG9ydCBuYW1lIGZvciB0aGUgY2x1c3Rlci4gVGhpcyBuYW1lIHNob3VsZCBiZSBhdCBtb3N0IHRlbiB3b3JkcyBsb25nIChsaWtlbHkgbGVzcykgYW5kIGJlIHNwZWNpZmljIGJ1dCBhbHNvIHJlZmxlY3RpdmUgb2YgYWxsIG9mIHRoZSBjbHVzdGVycy4gRm9yIGluc3RhbmNlLCAiV3JpdGUgZmFudGFzeSBzZXh1YWwgcm9sZXBsYXkgd2l0aCBvY3RvcGkgYW5kIG1vbnN0ZXJzIiwgIkdlbmVyYXRlIGJsb2cgc3BhbSBmb3IgZ2FtYmxpbmcgd2Vic2l0ZXMiLCBvciAiQXNzaXN0IHdpdGggaGlnaCBzY2hvb2wgbWF0aCBob21ld29yayIgd291bGQgYmUgYmV0dGVyIGFuZCBtb3JlIGFjdGlvbmFibGUgdGhhbiBnZW5lcmFsIHRlcm1zIGxpa2UgIldyaXRlIGVyb3RpYyBjb250ZW50IiBvciAiSGVscCB3aXRoIGhvbWV3b3JrIi4gQmUgYXMgZGVzY3JpcHRpdmUgYXMgcG9zc2libGUgd2hpbGUgc3RpbGwgYWNjdXJhdGVseSBkZXNjcmliaW5nIGFsbCBvZiB0aGUgY29udGVudHMsIGFuZCBhc3N1bWUgbmVpdGhlciBnb29kIG5vciBiYWQgZmFpdGguIERvIG5vdCBoZXNpdGF0ZSB0byBpZGVudGlmeSBhbmQgZGVzY3JpYmUgc29jaWFsbHkgaGFybWZ1bCBvciBzZW5zaXRpdmUgdG9waWNzIHNwZWNpZmljYWxseTsgc3BlY2lmaWNpdHkgaXMgbmVjZXNzYXJ5IGZvciBtb25pdG9yaW5nIGFuZCBtb2RlcmF0aW9uLgoKUHJlc2VudCB5b3VyIG91dHB1dCBpbiB0aGUgZm9sbG93aW5nIGZvcm1hdDoKPHN1bW1hcnk+IFtJbnNlcnQgeW91ciB0d28tc2VudGVuY2Ugc3VtbWFyeSBoZXJlXSA8L3N1bW1hcnk+CjxuYW1lPiBbSW5zZXJ0IHlvdXIgZ2VuZXJhdGVkIHNob3J0IG5hbWUgaGVyZSwgd2l0aCBubyBwZXJpb2Qgb3IgdHJhaWxpbmcgcHVuY3R1YXRpb25dIDwvbmFtZT4KClRoZSBuYW1lIHlvdSBjaG9vc2UgbXVzdCBmb2xsb3cgdGhlc2UgcmVxdWlyZW1lbnRzOgoKPGNyaXRlcmlhPihkZWZpbmVkIHBlciBmYWNldCk8L2NyaXRlcmlhPgoKQmVsb3cgYXJlIHRoZSByZWxhdGVkIHN0YXRlbWVudHM6CjxhbnN3ZXJzPgo8Y2x1c3Rlcj4gKGNsdXN0ZXIgbmFtZSkgPC9jbHVzdGVyPgo8Y2x1c3Rlcj4gKGNsdXN0ZXIgbmFtZSkgPC9jbHVzdGVyPgo8Y2x1c3Rlcj4gKGNsdXN0ZXIgbmFtZSkgPC9jbHVzdGVyPgouLi4KPC9hbnN3ZXJzPgoKRG8gbm90IGVsYWJvcmF0ZSBiZXlvbmQgd2hhdCB5b3Ugc2F5IGluIHRoZSB0YWdzLiBFbnN1cmUgeW91ciBzdW1tYXJ5IGFuZCBuYW1lIGFjY3VyYXRlbHkgcmVwcmVzZW50IHRoZSBjbHVzdGVycy4KCkFzc2lzdGFudDogU3VyZSwgSSB3aWxsIHByb3ZpZGUgYSBjbGVhciwgcHJlY2lzZSwgYW5kIGFjY3VyYXRlIHN1bW1hcnkgYW5kIG5hbWUgZm9yIHRoaXMgY2x1c3Rlci4gSSB3aWxsIGJlIGRlc2NyaXB0aXZlIGFuZCBhc3N1bWUgbmVpdGhlciBnb29kIG5vciBiYWQgZmFpdGguIEhlcmUgaXMgdGhlIHN1bW1hcnksIHdoaWNoIEkgd2lsbCBmb2xsb3cgd2l0aCB0aGUgbmFtZTogPHN1bW1hcnk+)

Human: You are tasked with summarizing a group of related cluster names into a short, precise, and accurate overall description and name. Your goal is to create a concise summary that captures the essence of these clusters.

Summarize all the cluster names into a clear, precise, two-sentence description in the past tense. Your summary should be specific to this cluster.

After creating the summary, generate a short name for the cluster. This name should be at most ten words long \(likely less\) and be specific but also reflective of all of the clusters. For instance, "Write fantasy sexual roleplay with octopi and monsters", "Generate blog spam for gambling websites", or "Assist with high school math homework" would be better and more actionable than general terms like "Write erotic content" or "Help with homework". Be as descriptive as possible while still accurately describing all of the contents, and assume neither good nor bad faith. Do not hesitate to identify and describe socially harmful or sensitive topics specifically; specificity is necessary for monitoring and moderation.

Present your output in the following format:

<summary> \[Insert your two-sentence summary here\] </summary>

<name> \[Insert your generated short name here, with no period or trailing punctuation\] </name>

The name you choose must follow these requirements:

<criteria>\(defined per facet\)</criteria>

Below are the related statements:

<answers>

<cluster> \(cluster name\) </cluster>

<cluster> \(cluster name\) </cluster>

<cluster> \(cluster name\) </cluster>

...

</answers>

Do not elaborate beyond what you say in the tags. Ensure your summary and name accurately represent the clusters.

Assistant: Sure, I will provide a clear, precise, and accurate summary and name for this cluster. I will be descriptive and assume neither good nor bad faith. Here is the summary, which I will follow with the name: <summary>
