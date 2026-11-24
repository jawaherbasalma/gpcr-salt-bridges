#!/usr/bin/env nextflow

/*
 * gpcr-salt-bridges -- the ionic lock pipeline as a Nextflow workflow.
 *
 * This runs exactly the same five steps as `pixi run reproduce`, in the same
 * order, against the same scripts. Nothing is duplicated: if a script changes,
 * both routes change with it. What Nextflow adds is an explicit dependency
 * graph, the ability to resume a half-finished run, and an execution report
 * recording how long each step took and how much memory it used.
 *
 *     pixi run nextflow
 *
 * To run the full structure set rather than the committed subset:
 *
 *     nextflow run workflow/main.nf --config config/structures_all.yaml -resume
 *
 * Why the steps run in the project directory
 * ------------------------------------------
 * Nextflow usually isolates each task in its own scratch directory and stages
 * input files into it. That is the right default, and it is deliberately not
 * what happens here.
 *
 * This project treats data/raw/ as a permanent store: structures are
 * downloaded once, their SHA256 checksums are recorded in
 * metadata/structures.csv, and the files are set read-only so that no later
 * step can modify them in place. Script 01 skips anything already present.
 * Staging those files into per-task scratch directories would mean either
 * re-downloading eleven structures from the RCSB on every run -- discourteous
 * to a public academic service, and slow -- or relaxing the read-only
 * guarantee that makes the provenance record meaningful.
 *
 * So each process changes into the project directory and runs its script
 * there, and the processes are chained by small completion tokens rather than
 * by staged data files. The ordering Nextflow enforces is therefore real; the
 * task-level sandboxing is not. That is the trade-off, made knowingly.
 *
 * The natural next step, if this grows
 * ------------------------------------
 * The slow step by far is arpeggio, at roughly a minute per structure, and
 * every structure is independent of every other. Splitting step 2 into one
 * task per structure would let Nextflow run them in parallel, which is the
 * point at which a workflow manager starts earning its keep rather than
 * merely documenting the order. That needs script 02 to accept a single PDB
 * ID, which it does not yet.
 */


process DOWNLOAD {
    tag 'structures'

    output:
    path 'download.done'

    script:
    """
    set -euo pipefail
    token=\$PWD/download.done
    cd ${projectDir}/..
    STRUCTURES_CONFIG='${params.config}' python scripts/01_download_structures.py
    touch \$token
    """
}


process SALT_BRIDGES {
    tag 'arpeggio'

    input:
    path done

    output:
    path 'salt_bridges.done'

    script:
    """
    set -euo pipefail
    token=\$PWD/salt_bridges.done
    cd ${projectDir}/..
    python scripts/02_extract_salt_bridges.py
    touch \$token
    """
}


process GENERIC_POSITIONS {
    tag 'gpcrdb'

    input:
    path done

    output:
    path 'positions.done'

    script:
    """
    set -euo pipefail
    token=\$PWD/positions.done
    cd ${projectDir}/..
    STRUCTURES_CONFIG='${params.config}' python scripts/03_map_generic_positions.py
    touch \$token
    """
}


process TEST_LOCK {
    tag 'statistics'

    input:
    path bridges_done
    path positions_done

    output:
    path 'test.done'

    script:
    """
    set -euo pipefail
    token=\$PWD/test.done
    cd ${projectDir}/..
    python scripts/04_test_lock.py
    touch \$token
    """
}


process FIGURES {
    tag 'figures'

    input:
    path done

    output:
    path 'figures.done'

    script:
    """
    set -euo pipefail
    token=\$PWD/figures.done
    cd ${projectDir}/..
    python scripts/figures.py
    touch \$token
    """
}


process REPORT {
    tag 'quarto'

    input:
    path done

    output:
    path 'report.done'

    script:
    """
    set -euo pipefail
    token=\$PWD/report.done
    cd ${projectDir}/..
    quarto render reports/report.qmd
    touch \$token
    """
}


workflow {

    // Steps 2 and 3 both depend only on step 1, and nothing on each other --
    // one reads the structures for contacts, the other for residue numbering.
    // Written this way, the graph records that honestly, and they would run
    // concurrently if the tasks were isolated.
    downloaded = DOWNLOAD()

    bridges   = SALT_BRIDGES(downloaded)
    positions = GENERIC_POSITIONS(downloaded)

    tested = TEST_LOCK(bridges, positions)
    drawn  = FIGURES(tested)
    REPORT(drawn)
}
