import logging
from typing import List, Optional
from datetime import datetime

import numpy as np

from . import helpers
from . import entrapment


logger = logging.getLogger(__name__)


def calculateProteinFDRs(proteinGroups, proteinScores):
    logger.info("Calculating protein group-level FDRs")
    numDecoys, numEntrapments, numTargets = 0, 0, 0
    proteinGroupInfoList = list()
    for proteinGroup, proteinScore in zip(proteinGroups, proteinScores):
        if proteinScore == -100.0:
            break

        if helpers.is_decoy(proteinGroup):
            numDecoys += 1
        else:
            numTargets += 1
            if entrapment.isEntrapment(proteinGroup):
                numEntrapments += 1
        reportedFdr = (numDecoys + 1) / (numTargets + 1)
        observedFdr = (numEntrapments + 1) / (numTargets + 1)

        skipForCounting = helpers.is_decoy(proteinGroup) or helpers.isObsolete(
            proteinGroup
        )
        proteinGroupInfoList.append((reportedFdr, observedFdr, skipForCounting))

    logger.info(f"#Protein groups: Target = {numTargets}; Decoys = {numDecoys}")
    if numEntrapments > 1:
        logger.info(
            f"    Entrapments = {numEntrapments}; Targets-Entrapments = {numTargets - numEntrapments}"
        )

    if len(proteinGroupInfoList) == 0:
        raise Exception(
            "No proteins with scores found, make sure that protein identifiers are consistent in the evidence and fasta files"
        )

    reportedFdrs, observedFdrs, skipForCounting = zip(*proteinGroupInfoList)
    reportedQvals, observedQvals = fdrsToQvals(reportedFdrs), fdrsToQvals(observedFdrs)
    logger.info(
        f"#Target protein groups at 1% decoy FDR: {countBelowThreshold(reportedQvals, 0.01, skipForCounting)}"
    )
    if numEntrapments > 1:
        logger.info(
            f"#Target protein groups at 1% entrapment FDR: {countBelowThreshold(observedFdrs, 0.01, skipForCounting)}"
        )
        num_decoys_at_entrapment_fdr = countBelowThreshold(observedFdrs, 0.01)
        reported_qvals_index = min(len(reportedQvals)-1, num_decoys_at_entrapment_fdr)
        logger.info(
            f"Decoy FDR at 1% entrapment FDR: {'%.2g' % (reportedQvals[reported_qvals_index])}"
        )
        num_entrapments_at_decoy_fdr = countBelowThreshold(reportedQvals, 0.01)
        observed_fdrs_index = min(len(observedFdrs)-1, num_entrapments_at_decoy_fdr)
        logger.info(
            f"Entrapment FDR at 1% decoy FDR: {'%.2g' % (observedFdrs[observed_fdrs_index])}"
        )

        # printReportedAndEntrapmentFDRs(reportedQvals, observedQvals)

    return reportedQvals, observedQvals


def calculatePeptideFDRs(peptideScores, scoreType):
    decoyScores, entrapmentScores, poolScores = list(), list(), list()
    fdrs = list()
    reportedFdr, observedFdr = 0, 0
    sumPEP = np.nextafter(0, 1)
    for score, _, proteins in peptideScores:
        if helpers.is_contaminant(proteins):
            continue

        if helpers.is_decoy(proteins):
            decoyScores.append(score)
        else:
            if entrapment.isEntrapment(proteins):
                entrapmentScores.append(score)
            else:
                poolScores.append(score)
            if "PEPavg" in scoreType.long_description():
                sumPEP += score
                reportedFdr = sumPEP / (len(entrapmentScores) + len(poolScores))
            else:
                reportedFdr = (len(decoyScores) + 1) / (
                    len(poolScores) + len(entrapmentScores) + 1
                )
            observedFdr = (len(entrapmentScores) + 1) / (
                len(poolScores) + len(entrapmentScores) + 1
            )
            fdrs.append((reportedFdr, observedFdr))
        # print(score, reportedFdr, observedFdr, helpers.isDecoy(proteins), entrapment.isEntrapment(proteins), sep='\t')

    numDecoys = len(decoyScores)
    numEntrapments = len(entrapmentScores)
    numTargets = len(poolScores) + len(entrapmentScores)
    logger.info(f"#Peptides: Target = {numTargets}; Decoys = {numDecoys}")
    if numEntrapments > 1:
        logger.info(
            f"    Entrapments = {numEntrapments}; Targets-Entrapments = {numTargets - numEntrapments}"
        )

    reportedFdrs, observedFdrs = zip(*fdrs)
    reportedQvals, observedQvals = fdrsToQvals(reportedFdrs), fdrsToQvals(observedFdrs)
    logger.info(
        f"#Target peptides at 1% decoy FDR: {countBelowThreshold(reportedQvals, 0.01)}"
    )
    if numEntrapments > 1:
        logger.info(
            f"#Target peptides at 1% entrapment FDR: {countBelowThreshold(observedFdrs, 0.01)}"
        )
        num_decoys_at_entrapment_fdr = countBelowThreshold(observedFdrs, 0.01)
        reported_qvals_index = min(len(reportedQvals)-1, num_decoys_at_entrapment_fdr)
        logger.info(
            f"Decoy FDR at 1% entrapment FDR: {'%.2g' % (reportedQvals[reported_qvals_index])}"
        )
        num_entrapments_at_decoy_fdr = countBelowThreshold(reportedQvals, 0.01)
        observed_fdrs_index = min(len(observedFdrs)-1, num_entrapments_at_decoy_fdr)
        logger.info(
            f"Entrapment FDR at 1% decoy FDR: {'%.2g' % (observedFdrs[observed_fdrs_index])}"
        )

    return reportedQvals, observedQvals


def printReportedAndEntrapmentFDRs(reportedQvals, observedQvals):
    import csv

    writer = csv.writer(
        open(
            f'protein_fdr_calibration_{datetime.now().strftime("%d%m%Y_%H%M%S")}.txt',
            "w",
        ),
        delimiter="\t",
    )
    for reportedQval, observedQval in zip(reportedQvals, observedQvals):
        writer.writerow([reportedQval, observedQval])


def fdrsToQvals(fdrs: List[float]) -> np.array:
    """
    Makes a list of FDRs monotonically increasing (sometimes referred to as q-values after monotonization)
    """
    return np.minimum.accumulate(fdrs[::-1])[::-1]


def countBelowThreshold(
    qvals: List[float],
    qvalThreshold: float,
    skipForCounting: Optional[List[bool]] = None,
):
    """
    Counts number of q-values below a threshold, if skipForCounting are provided, only the targets are counted
    """
    if skipForCounting is None:
        return len([1 for x in qvals if x < qvalThreshold])
    else:
        return len(
            [
                1
                for x, skip in zip(qvals, skipForCounting)
                if x < qvalThreshold and not skip
            ]
        )


def calc_post_err_prob_cutoff(postErrProbs, psmQvalCutoff):
    postErrProbCutoff = 1.0
    sumPEP = 0.0
    numPSMs = 0
    for postErrProb in sorted(postErrProbs):
        if not np.isfinite(postErrProb):
            continue
        sumPEP += postErrProb
        numPSMs += 1
        if sumPEP / numPSMs > psmQvalCutoff:
            postErrProbCutoff = postErrProb
            break
    return postErrProbCutoff
