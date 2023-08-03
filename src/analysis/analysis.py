### Run analysis

columns_chembl = ["LoF_protect", "GoF_protect"]
columns_dataset = ["LoF_protect", "GoF_protect", "LoF_risk", "GoF_risk", "evidenceDif"]
columns = ["GoF_risk", "LoF_protect", "LoF_risk", "GoF_protect"]
terms = ["noEvaluable", "bivalent_risk", "null", "dispar"]

sincgc = [
    "gene_burden",
    "intogen",
    "eva",
    "eva_somatic",
    "ot_genetics_portal",
    "impc",
    "orphanet",
    "gene2phenotype",
]

germline = [
    "gene_burden",
    "eva",
    "ot_genetics_portal",
    "impc",
    "orphanet",
    "gene2phenotype",
]

somatic = ["intogen", "cancer_gene_census", "eva_somatic"]

datasource_list = [
    "gene_burden",
    "intogen",
    "cancer_gene_census",
    "eva",
    "eva_somatic",
    "ot_genetics_portal",
    "impc",
    "orphanet",
    "gene2phenotype",
    "chembl",
    "WOcgc",
    "somatic",
    "germline",
]

genEvidDataset = (
    prueba_assessment.filter(F.col("datasourceId") != "chembl")  #### checked 31.05.2023
    .groupBy("targetId", "diseaseId")
    .agg(F.count("targetId").alias("Nr_evidences"))
    .select("targetId", "diseaseId", "Nr_evidences")
    .withColumn("geneticEvidence", F.lit("hasGeneticEvidence"))
)

coherency_toAssess_others_datasource = (  #### checked 31.05.2023
    prueba_assessment.filter(
        (F.col("Assessment").isin(terms) == False) & (F.col("datasourceId") != "chembl")
    )
    .groupBy("targetId", "diseaseId")
    .agg(F.collect_set("datasourceId").alias("datasourceIds"))
)

childrens = diseases.withColumn(
    "propagatedTrait",
    F.explode_outer(F.concat(F.array(F.col("id")), F.col("descendants"))),
).select("id", "propagatedTrait")

wByDisease = Window.partitionBy("diseaseId")  #### checked 31.05.2023
diseaseTA = (
    diseases.withColumn("taId", F.explode("therapeuticAreas"))
    .select(F.col("id").alias("diseaseId"), "taId")
    .join(taDf, on="taId", how="left")
    .withColumn("minRank", F.min("taRank").over(wByDisease))
    .filter(F.col("taRank") == F.col("minRank"))
    .drop("taRank", "minRank")
)

childrensTA = childrens.join(
    diseaseTA.select("diseaseId", "taLabelSimple"),
    F.col("propagatedTrait") == diseaseTA.diseaseId,
    "left",
)

taDf = spark.createDataFrame(
    data=[
        ("MONDO_0045024", "cell proliferation disorder", "Oncology"),
        ("EFO_0005741", "infectious disease", "Other"),
        ("OTAR_0000014", "pregnancy or perinatal disease", "Other"),
        ("EFO_0005932", "animal disease", "Other"),
        ("MONDO_0024458", "disease of visual system", "Other"),
        ("EFO_0000319", "cardiovascular disease", "Other"),
        ("EFO_0009605", "pancreas disease", "Other"),
        ("EFO_0010282", "gastrointestinal disease", "Other"),
        ("OTAR_0000017", "reproductive system or breast disease", "Other"),
        ("EFO_0010285", "integumentary system disease", "Other"),
        ("EFO_0001379", "endocrine system disease", "Other"),
        ("OTAR_0000010", "respiratory or thoracic disease", "Other"),
        ("EFO_0009690", "urinary system disease", "Other"),
        ("OTAR_0000006", "musculoskeletal or connective tissue disease", "Other"),
        ("MONDO_0021205", "disease of ear", "Other"),
        ("EFO_0000540", "immune system disease", "Other"),
        ("EFO_0005803", "hematologic disease", "Other"),
        ("EFO_0000618", "nervous system disease", "Other"),
        ("MONDO_0002025", "psychiatric disorder", "Other"),
        ("MONDO_0024297", "nutritional or metabolic disease", "Other"),
        ("OTAR_0000018", "genetic, familial or congenital disease", "Other"),
        ("OTAR_0000009", "injury, poisoning or other complication", "Other"),
        ("EFO_0000651", "phenotype", "Other"),
        ("EFO_0001444", "measurement", "Other"),
        ("GO_0008150", "biological process", "Other"),
    ],
    schema=StructType(
        [
            StructField("taId", StringType(), True),
            StructField("taLabel", StringType(), True),
            StructField("taLabelSimple", StringType(), True),
        ]
    ),
).withColumn("taRank", F.monotonically_increasing_id())


chemblPhase = chembl.groupBy("targetId", "diseaseId").agg(
    F.max("clinicalPhase").alias("maxClinPhase")
)

chembl_forGenEvid = (
    prueba_assessment.filter(F.col("datasourceId") == "chembl")  #### checked 31.05.2023
    .groupBy("targetId", "diseaseId", "datasourceId")
    .pivot("homogenized")
    .agg(F.count("targetId"))
    .select(
        F.col("targetId").alias("targetIdC"),
        F.col("datasourceId").alias("datasourceIdC"),
        F.col("diseaseId").alias("diseaseIdC"),
        *(F.col(c).cast("int").alias(c) for c in columns_chembl),
    )
    .withColumn(
        "coherency_chembl",  #### checked 31.05.2023
        F.when(
            F.col("GoF_protect").isNotNull(),
            F.when(F.col("LoF_protect").isNotNull(), F.lit("dispar")).when(
                F.col("LoF_protect").isNull(), F.lit("coherent")
            ),
        )
        .when(
            F.col("LoF_protect").isNotNull(),
            F.when(F.col("GoF_protect").isNotNull(), F.lit("dispar")).when(
                F.col("GoF_protect").isNull(), F.lit("coherent")
            ),
        )
        .when(
            ((F.col("GoF_protect").isNull()) & (F.col("LoF_protect").isNull())),
            F.lit("noData"),
        ),
    )
    .join(
        chemblPhase,
        (chemblPhase.targetId == F.col("targetIdC"))
        & (chemblPhase.diseaseId == F.col("diseaseIdC")),
        "left",
    )
    .drop("targetId", "diseaseId")
    .join(
        genEvidDataset,
        (genEvidDataset.targetId == F.col("targetIdC"))
        & (genEvidDataset.diseaseId == F.col("diseaseIdC")),
        "left",
    )
    .drop("targetId", "diseaseId")
    .join(
        diseaseTA.select("diseaseId", "taLabelSimple"),
        F.col("diseaseIdC") == diseaseTA.diseaseId,
        "left",
    )
    .drop("diseaseId")
    .withColumn(
        "geneticEvidenceAssessment",
        F.when(F.col("geneticEvidence").isNotNull(), F.lit("hasGenetics")).otherwise(
            F.lit("noGenetics")
        ),
    )
    .withColumnRenamed("LoF_protect", "LoF_protectC")
    .withColumnRenamed("GoF_protect", "GoF_protectC")
    .withColumnRenamed("maxClinPhase", "maxClinPhase_ext")
)


def builder_original_genEvid(  #### checked 31.05.2023
    df, columns_dataset, chembl_forGenEvid
):  #### checked 31.05.2023
    df = (
        df.select(
            F.col("targetId").alias("targetIdO"),
            F.col("diseaseId").alias("diseaseIdO"),
            *(F.col(c).cast("int").alias(c) for c in columns_dataset),
        )
        .groupBy("targetIdO", "diseaseIdO")
        .agg(
            F.sum("GoF_protect").alias("GoF_protect"),
            F.sum("LoF_protect").alias("LoF_protect"),
            F.sum("LoF_risk").alias("LoF_risk"),
            F.sum("GoF_risk").alias("GoF_risk"),
            F.sum("evidenceDif").alias("evidenceDif"),
        )
        .withColumn(
            "coherency_others",  #### checked 31.05.2023
            F.when(
                (
                    F.col("GoF_risk").isNotNull()
                    | F.col("LoF_risk").isNotNull()
                    | F.col("LoF_protect").isNotNull()
                    | F.col("GoF_protect").isNotNull()
                ),
                F.when(
                    ((F.col("GoF_risk").isNotNull()) & (F.col("LoF_risk").isNotNull())),
                    F.lit("dispar"),
                )
                .when(
                    (
                        (F.col("LoF_protect").isNotNull())
                        & (F.col("LoF_risk").isNotNull())
                    ),
                    F.lit("dispar"),
                )
                .when(
                    (
                        (F.col("GoF_protect").isNotNull())
                        & (F.col("GoF_risk").isNotNull())
                    ),
                    F.lit("dispar"),
                )
                .when(
                    (
                        (F.col("GoF_protect").isNotNull())
                        & (F.col("LoF_protect").isNotNull())
                    ),
                    F.lit("dispar"),
                )
                .otherwise(F.lit("coherent")),
            ).otherwise(F.lit("evidenceDif")),
        )
        .join(  ### checked 01.06.2023
            chembl_forGenEvid,
            (F.col("targetIdO") == chembl_forGenEvid.targetIdC)
            & (F.col("diseaseIdO") == chembl_forGenEvid.diseaseIdC),
            "right",
        )
        .withColumn(  #### checked 01.06.2023
            "coherency_inter",
            F.when(
                F.col("targetIdO").isNotNull(),
                F.when(
                    (
                        (F.col("coherency_others") == "coherent")
                        & (F.col("coherency_chembl") == "coherent")
                    ),
                    F.when(
                        (F.col("GoF_protectC").isNotNull())
                        & (
                            (F.col("LoF_protect").isNotNull())
                            | (F.col("GoF_risk").isNotNull())
                        ),
                        F.lit("noDoE/noGenetics"),
                    )
                    .when(
                        (F.col("LoF_protectC").isNotNull())
                        & (
                            (F.col("GoF_protect").isNotNull())
                            | (F.col("LoF_risk").isNotNull())
                        ),
                        F.lit("noDoE/noGenetics"),
                    )
                    .otherwise(F.lit("coherentDoE")),
                ).otherwise(F.lit("noDoE/noGenetics")),
            )
            .when(
                F.col("targetIdO").isNull(),
                F.lit("noDoE/noGenetics"),
            )
            .otherwise(F.lit("noDoE/noGenetics")),
        )
        .withColumn(  #### checked 01.06.2023
            "coherency_inter_oneCell",
            F.when(
                F.col("targetIdO").isNotNull(),
                F.when(
                    (
                        (F.col("coherency_others") == "coherent")
                        & (F.col("coherency_chembl") == "coherent")
                    ),
                    F.when(
                        (F.col("GoF_protectC").isNotNull())
                        & (F.col("GoF_protect").isNotNull()),
                        F.lit("coherentDoE"),
                    )
                    .when(
                        (F.col("LoF_protectC").isNotNull())
                        & (F.col("LoF_protect").isNotNull()),
                        F.lit("coherentDoE"),
                    )
                    .otherwise(F.lit("noDoE/noGenetics")),
                ).otherwise(F.lit("noDoE/noGenetics")),
            )
            .when(
                F.col("targetIdO").isNull(),
                F.lit("noDoE/noGenetics"),
            )
            .otherwise(F.lit("noDoE/noGenetics")),
        )
        #### two columns for the direction of effect (over the ones with DoE)
        .withColumn(  #### checked 01.06.2023
            "coherencyOnlyDOE_inter",
            F.when(
                F.col("targetIdO").isNotNull(),
                F.when(
                    (
                        (F.col("coherency_others") == "coherent")
                        & (F.col("coherency_chembl") == "coherent")
                    ),
                    F.when(
                        (F.col("GoF_protectC").isNotNull())
                        & (
                            (F.col("LoF_protect").isNotNull())
                            | (F.col("GoF_risk").isNotNull())
                        ),
                        F.lit("dispar"),
                    )
                    .when(
                        (F.col("LoF_protectC").isNotNull())
                        & (
                            (F.col("GoF_protect").isNotNull())
                            | (F.col("LoF_risk").isNotNull())
                        ),
                        F.lit("dispar"),
                    )
                    .otherwise(F.lit("coherent")),
                ).otherwise(F.lit(None)),
            )
            .when(
                F.col("targetIdO").isNull(),
                F.lit(None),
            )
            .otherwise(F.lit(None)),
        )
        .withColumn(  #### checked 01.06.2023
            "coherencyOnlyDOE_inter_oneCell",
            F.when(
                F.col("targetIdO").isNotNull(),
                F.when(
                    (
                        (F.col("coherency_others") == "coherent")
                        & (F.col("coherency_chembl") == "coherent")
                    ),
                    F.when(
                        (F.col("GoF_protectC").isNotNull())
                        & (F.col("GoF_protect").isNotNull()),
                        F.lit("coherent"),
                    )
                    .when(
                        (F.col("LoF_protectC").isNotNull())
                        & (F.col("LoF_protect").isNotNull()),
                        F.lit("coherent"),
                    )
                    .otherwise(F.lit("dispar")),
                ).otherwise(F.lit(None)),
            )
            .when(
                F.col("targetIdO").isNull(),
                F.lit(None),
            )
            .otherwise(F.lit(None)),
        )
        ### genetic Evidence from every datasource & indication
        .withColumn(
            "hasGeneticsFromDS",
            F.when(
                F.col("targetIdO").isNotNull(), F.lit("GeneticEvidenceFromDS")
            ).otherwise(F.lit("noGeneticEvidence")),
        )
        .withColumn(
            "oncologyHasGeneticsFromDS",
            F.when(
                (F.col("targetIdO").isNotNull())
                & (F.col("taLabelSimple") == "Oncology"),
                F.lit("hasGeneticEvidence"),
            )
            .when(
                (F.col("targetIdO").isNull()) & (F.col("taLabelSimple") == "Oncology"),
                F.lit("noGeneticEvidence"),
            )
            .otherwise(F.lit(None)),
        )
        .withColumn(
            "otherHasGeneticsFromDS",
            F.when(
                (F.col("targetIdO").isNotNull()) & (F.col("taLabelSimple") == "Other"),
                F.lit("hasGeneticEvidence"),
            )
            .when(
                (F.col("targetIdO").isNull()) & (F.col("taLabelSimple") == "Other"),
                F.lit("noGeneticEvidence"),
            )
            .otherwise(F.lit(None)),
        )
        ##oncology & genetics
        .withColumn(  #### checked 01.06.2023
            "oncologyInd_genetics",
            F.when(
                F.col("taLabelSimple") == "Oncology", F.col("geneticEvidenceAssessment")
            ).otherwise(F.lit(None)),
        )
        ## oncology & coherency_inter
        .withColumn(  #### checked 01.06.2023
            "oncologyInd_coherency_inter",
            F.when(
                F.col("taLabelSimple") == "Oncology", F.col("coherency_inter")
            ).otherwise(F.lit(None)),
        )
        ## oncology & coherency_interOneCell
        .withColumn(  #### checked 01.06.2023
            "oncologyInd_coherency_interOneCell",
            F.when(
                F.col("taLabelSimple") == "Oncology",
                F.col("coherency_inter_oneCell"),
            ).otherwise(F.lit(None)),
        )
        ## other & genetics
        .withColumn(  #### checked 01.06.2023
            "otherInd_genetics",
            F.when(
                F.col("taLabelSimple") == "Other", F.col("geneticEvidenceAssessment")
            ).otherwise(F.lit(None)),
        )
        ## other & coherency_inter
        .withColumn(  #### checked 01.06.2023
            "otherInd_coherency_inter",
            F.when(
                F.col("taLabelSimple") == "Other", F.col("coherency_inter")
            ).otherwise(F.lit(None)),
        )
        ## other & coherency_interOneCell
        .withColumn(  #### checked 01.06.2023
            "otherInd_coherency_interOneCell",
            F.when(
                F.col("taLabelSimple") == "Other", F.col("coherency_inter_oneCell")
            ).otherwise(F.lit(None)),
        )
        #### oncology only DOE coherencyInter
        .withColumn(  #### checked 02.06.2023
            "oncologyIndOnlyDoe_coherency_inter",
            F.when(
                F.col("taLabelSimple") == "Oncology", F.col("coherencyOnlyDOE_inter")
            ).otherwise(F.lit(None)),
        )
        #### oncology only DOE coherencyOneCell
        .withColumn(  #### checked  02.06.2023
            "oncologyIndOnlyDoe_coherency_interOneCell",
            F.when(
                F.col("taLabelSimple") == "Oncology",
                F.col("coherencyOnlyDOE_inter_oneCell"),
            ).otherwise(F.lit(None)),
        )
        #### other only DOE coherencyInter
        .withColumn(  #### checked  02.06.2023
            "otherIndOnlyDoe_coherency_inter",
            F.when(
                F.col("taLabelSimple") == "Other", F.col("coherencyOnlyDOE_inter")
            ).otherwise(F.lit(None)),
        )
        #### other only DOE coherencyOneCell
        .withColumn(  #### checked  02.06.2023
            "otherIndOnlyDoe_coherency_interOneCell",
            F.when(
                F.col("taLabelSimple") == "Other",
                F.col("coherencyOnlyDOE_inter_oneCell"),
            ).otherwise(F.lit(None)),
        )
        ## phase clin trials
        .withColumn(  #### checked 01.06.2023
            "pro_phase4",
            F.when(F.col("maxClinPhase_ext") == 4, F.lit("PhaseIV")).otherwise(
                F.lit("dif")
            ),
        )
        .withColumn(  #### checked 01.06.2023
            "pro_phase3",
            F.when(F.col("maxClinPhase_ext") == 3, F.lit("PhaseIII")).otherwise(
                F.lit("dif")
            ),
        )
        .withColumn(  #### checked 01.06.2023
            "pro_phase2",
            F.when(F.col("maxClinPhase_ext") == 2, F.lit("PhaseII")).otherwise(
                F.lit("dif")
            ),
        )
        .withColumn(  #### checked 01.06.2023
            "pro_phase1",
            F.when(F.col("maxClinPhase_ext") == 1, F.lit("PhaseI")).otherwise(
                F.lit("dif")
            ),
        )
        .withColumn(  #### checked 01.06.2023
            "pro_phase0",
            F.when(F.col("maxClinPhase_ext") == 0.5, F.lit("Phase0")).otherwise(
                F.lit("dif")
            ),
        )
        .withColumn(  #### checked 01.06.2023
            "pro_phase3H",
            F.when(F.col("maxClinPhase_ext") >= 3, F.lit("PhaseIII+")).otherwise(
                F.lit("dif")
            ),
        )
        .withColumn(  #### checked 01.06.2023
            "pro_phase2H",
            F.when(F.col("maxClinPhase_ext") >= 2, F.lit("PhaseII+")).otherwise(
                F.lit("dif")
            ),
        )
        .withColumn(  #### checked 01.06.2023
            "pro_phase1H",
            F.when(
                F.col("maxClinPhase_ext") >= 1, F.lit("PhaseI+")
            ).otherwise(  ## fixed II+ 03.08.23
                F.lit("dif")
            ),
        )
    ).persist()
    return df


## 1nd dictionary
dfs_dict_genEvi_more = {}  ### checked and changed on 01.06.2023

assessment = prueba_assessment.withColumn(
    "unified",
    F.when(F.col("Assessment").isin(terms), F.lit("evidenceDif"))
    .when(F.col("Assessment").isNull(), F.lit("evidenceDif"))
    .when(F.col("Assessment") == "KO_risk", F.lit("LoF_risk"))
    .otherwise(F.col("homogenized")),
)

for value in datasource_list:

    if value == "chembl":
        filtered_df = (
            assessment.filter(F.col("datasourceId") != value)
            .groupBy("targetId", "diseaseId")
            .pivot("unified")
            .agg(F.count("targetId"))
        )

        key_name = f"df_!{value}"
        dfs_dict_genEvi_more[key_name] = filtered_df

    elif value == "WOcgc":
        filtered_df = (
            assessment.filter(F.col("datasourceId").isin(sincgc))
            .groupBy("targetId", "diseaseId")
            .pivot("unified")
            .agg(F.count("targetId"))
        )

        key_name = f"df_{value}"
        dfs_dict_genEvi_more[key_name] = filtered_df

    elif value == "germline":
        filtered_df = (
            assessment.filter(F.col("datasourceId").isin(germline))
            .groupBy("targetId", "diseaseId")
            .pivot("unified")
            .agg(F.count("targetId"))
        )

        key_name = f"df_{value}"
        dfs_dict_genEvi_more[key_name] = filtered_df

    elif value == "somatic":
        filtered_df = (
            assessment.filter(F.col("datasourceId").isin(somatic))
            .groupBy("targetId", "diseaseId")
            .pivot("unified")
            .agg(F.count("targetId"))
        )

        key_name = f"df_{value}"
        dfs_dict_genEvi_more[key_name] = filtered_df

    else:

        filtered_df = (
            assessment.filter(F.col("datasourceId") == value)
            .groupBy("targetId", "diseaseId")
            .pivot("unified")
            .agg(F.count("targetId"))
        )

        key_name = f"df_{value}"
        dfs_dict_genEvi_more[key_name] = filtered_df


#### original target-trait-drug pairs #### checked 01.06.2023
for key, df in dfs_dict_genEvi_more.items():

    for col in columns_dataset:
        if col not in df.columns:
            df = df.withColumn(col, F.lit(None))

    df = builder_original_genEvid(df, columns_dataset, chembl_forGenEvid)

    dfs_dict_genEvi_more[key] = df

### 26.05.2023
### analysis Original
comparisons_dataSource_original = spark.createDataFrame(
    data=[
        ("coherency_inter", "byDatatype"),
        ("coherency_inter_oneCell", "byDatatype"),
        # ("geneticEvidenceAssessment","byDataType"), ### same data as hasGeneticsFromDS
        # ("oncologyInd_genetics", "byDatatype"),
        ("oncologyInd_coherency_inter", "byDatatype"),
        ("oncologyInd_coherency_interOneCell", "byDataType"),
        # ("otherInd_genetics", "byDatatype"),
        ("otherInd_coherency_inter", "byDatatype"),
        ("otherInd_coherency_interOneCell", "byDataType"),
        ("coherencyOnlyDOE_inter", "byDataType"),
        ("coherencyOnlyDOE_inter_oneCell", "byDataType"),
        ("oncologyIndOnlyDoe_coherency_inter", "byDataType"),
        ("oncologyIndOnlyDoe_coherency_interOneCell", "byDataType"),
        ("otherIndOnlyDoe_coherency_inter", "byDataType"),
        ("otherIndOnlyDoe_coherency_interOneCell", "byDataType"),
        ("hasGeneticsFromDS", "byDataType"),
        ("oncologyHasGeneticsFromDS", "byDataType"),
        ("otherHasGeneticsFromDS", "byDataType"),
    ],
    schema=StructType(
        [
            StructField("comparison", StringType(), True),
            StructField("comparisonType", StringType(), True),
        ]
    ),
)

predictions_dataSource_original = spark.createDataFrame(
    data=[
        ("pro_phase4", "clinical"),
        ("pro_phase3", "clinical"),
        ("pro_phase2", "clinical"),
        ("pro_phase1", "clinical"),
        ("pro_phase0", "clinical"),
        ("pro_phase3H", "clinical"),
        ("pro_phase2H", "clinical"),
        ("pro_phase1H", "clinical"),
    ]
)
### function to write the data properly as a new folder (run before the above chunk)


def aggregations_original(
    df,
    data,
    listado_original,
    comparisonColumn,
    comparisonType,
    predictionColumn,
    predictionType,
):
    wComparison = Window.partitionBy(comparisonColumn)
    wPrediction = Window.partitionBy(predictionColumn)
    wPredictionComparison = Window.partitionBy(comparisonColumn, predictionColumn)
    ### original to join :
    original = df
    uniqIds = original.select("targetIdC", "diseaseIdC").distinct().count()
    ### los datos que dan evidencia genetica
    out = (
        original.withColumn("comparisonType", F.lit(comparisonType))
        .withColumn("predictionType", F.lit(predictionType))
        .withColumn("total", F.lit(uniqIds))
        .withColumn("a", F.count("targetIdC").over(wPredictionComparison))
        .withColumn(
            "predictionTotal",
            F.count("targetIdC").over(wPrediction),
        )
        .withColumn(
            "comparisonTotal",
            F.count("targetIdC").over(wComparison),
        )
        .select(
            F.col(predictionColumn).alias("prediction"),
            F.col(comparisonColumn).alias("comparison"),
            "comparisonType",
            "predictionType",
            "a",
            "predictionTotal",
            "comparisonTotal",
            "total",
        )
        .filter(F.col("prediction").isNotNull())
        .filter(F.col("comparison").isNotNull())
        .distinct()
    )

    today_date = str(date.today())
    out.write.mode("overwrite").parquet(
        today_date
        + "_"
        + "original/"
        + "_"
        + data
        + "_"
        + comparisonColumn
        + "_"
        + predictionColumn
        + ".parquet"
    )
    listado_original.append(
        str(
            today_date
            + "_"
            + "original/"
            + "_"
            + data
            + "_"
            + comparisonColumn
            + "_"
            + predictionColumn
            + ".parquet"
        )
    )
    ###print("datasources/" + "_" + data + "_" + comparisonColumn + "_" + predictionColumn + ".parquet")


aggSetups_original = comparisons_dataSource_original.join(
    predictions_dataSource_original, how="full"
).collect()
listado_original = []
for key, df in listado_original.items():
    for row in aggSetups_original:
        aggregations_original(df, key, listado_original, *row)
