PICKEDGROUPFDR_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include $(PICKEDGROUPFDR_DIR)MakefileShared

dependencies:
	git config --global credential.helper cache

registry:
	docker login gitlab.lrz.de:5005
	docker build -t gitlab.lrz.de:5005/proteomics/picked_group_fdr .
	docker push gitlab.lrz.de:5005/proteomics/picked_group_fdr

jump: 
	$(DOCKER_CMD) \
		$(IMAGE) bash

test:
	python3 -m pytest --cov=picked_group_fdr --cov-report html --cov-report term tests/unit_tests

performance:
	python3 -m pytest -s tests/performance_tests

line_profiler_lfq:
	kernprof -lv tests/performance_tests/test_lfq.py

memory_profile:
	mprof run --include-children --backend psutil_pss python3 -u tests/performance_tests/test_lfq.py | ts '[%H:%M:%.S]'

# --no-cache
build: dependencies
	docker build -f Dockerfile -t $(IMAGE) . || (exit 1)


prepayload_setup_create_folder: rm_err_file
	$(DOCKER_CMD) \
		$(IMAGE) mkdir -p -m 777 $(OUT_DIR_LOCAL)/percolator || (echo "1" > $(DATA)err.out; exit 1)

tab: prepayload_setup_create_folder
ifeq ($(PROSIT_FLAG),)
	$(DOCKER_CMD) \
		$(IMAGE) bash -c "python3 -u -m picked_group_fdr.pipeline.andromeda2pin $(MQ_EVIDENCE_FILE) --outputTab $(OUT_DIR_LOCAL)/percolator/andromeda.tab --databases \"$(LOCAL_DIR)/$(FASTA_FILE)\" $(DIGEST_PARAMS) > $(OUT_DIR_LOCAL)/percolator/andromeda2pin.log" || (echo "2" > $(DATA)err.out; exit 2)
endif

# calculating PEPs with Triqler within mokapot is very slow with OMP multithreading, so set this to 1 for now
percolator: tab
ifeq ($(PROSIT_FLAG),)
	$(DOCKER_CMD) \
		$(IMAGE) bash -c "OMP_NUM_THREADS=1 python3 -u -m picked_group_fdr.pipeline.run_mokapot $(PERC_TEST_FDR) $(PERC_TRAIN_FDR) $(OUT_DIR_LOCAL)/percolator" || (echo "3" > $(DATA)err.out; exit 3)
endif

update_evidence: percolator
	$(DOCKER_CMD) \
		$(IMAGE) bash -c "python3 -u -m picked_group_fdr.pipeline.update_evidence_from_pout --mq_evidence $(MQ_EVIDENCE_FILE) --perc_results $(PERC_RESULT_FILES) --mq_evidence_out $(OUT_DIR_LOCAL)/percolator/evidence.txt $(PROSIT_FLAG) > $(OUT_DIR_LOCAL)/percolator/update_evidence.log" || (echo "4" > $(DATA)err.out; exit 4)

picked_fdr: update_evidence
	$(DOCKER_CMD) \
		$(IMAGE) bash -c "OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 python3 -u -m picked_group_fdr --mq_evidence $(OUT_DIR_LOCAL)/percolator/evidence.txt --protein_groups_out $(OUT_DIR_LOCAL)/percolator/proteinGroups.txt --do_quant --fasta \"$(LOCAL_DIR)/$(FASTA_FILE)\" --num_threads $(NUM_THREADS) $(DIGEST_PARAMS) $(PICKED_GROUP_FDR_EXTRA_PARAMS) > $(OUT_DIR_LOCAL)/percolator/proteinGroups.log" || (echo "5" > $(DATA)err.out; exit 5)

picked_fdr_gene_level: picked_fdr
	$(DOCKER_CMD) \
		$(IMAGE) bash -c "OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 python3 -u -m picked_group_fdr --mq_evidence $(OUT_DIR_LOCAL)/percolator/evidence.txt --protein_groups_out $(OUT_DIR_LOCAL)/percolator/geneGroups.txt --do_quant --gene_level --fasta \"$(LOCAL_DIR)/$(FASTA_FILE)\" --num_threads $(NUM_THREADS) $(DIGEST_PARAMS) $(PICKED_GROUP_FDR_EXTRA_PARAMS) --suppress_missing_peptide_warning > $(OUT_DIR_LOCAL)/percolator/geneGroups.log" || (echo "5" > $(DATA)err.out; exit 5)

# by popular demand: also produce an evidence.txt and proteinGroups.txt filtered at 1% FDR
filter_results: picked_fdr_gene_level
	$(DOCKER_CMD) \
		$(IMAGE) bash -c "OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 python3 -u -m picked_group_fdr.pipeline.filter_fdr_maxquant --mq_msms $(OUT_DIR_LOCAL)/percolator/evidence.txt --mq_msms_out $(OUT_DIR_LOCAL)/percolator/evidence_fdr0.01.txt --mq_protein_groups $(OUT_DIR_LOCAL)/percolator/proteinGroups.txt --mq_protein_groups_out $(OUT_DIR_LOCAL)/percolator/proteinGroups_fdr0.01.txt --fdr_cutoff 0.01 --psm_level_fdr > $(OUT_DIR_LOCAL)/percolator/filter_results.log" || (echo "6" > $(DATA)err.out; exit 6)

filter_results_gene_level: filter_results
	$(DOCKER_CMD) \
		$(IMAGE) bash -c "python3 -u -m picked_group_fdr.pipeline.filter_fdr_maxquant --mq_protein_groups $(OUT_DIR_LOCAL)/percolator/geneGroups.txt --mq_protein_groups_out $(OUT_DIR_LOCAL)/percolator/geneGroups_fdr0.01.txt --fdr_cutoff 0.01 --psm_level_fdr > $(OUT_DIR_LOCAL)/percolator/filter_results_gene_level.log" || (echo "6" > $(DATA)err.out; exit 6)

compress: filter_results_gene_level
	zip -j -r -9 "$(OUT_DIR)/results.zip" "$(OUT_DIR)/percolator/" || (echo "" > $(DATA)err.out; exit 7)

all: compress


