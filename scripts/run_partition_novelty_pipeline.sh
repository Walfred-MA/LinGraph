#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage:
  bash minsetref/run_partition_novelty_pipeline.sh [OPTIONS]
  bash minsetref/run_partition_novelty_pipeline.sh [OPTIONS] [COHORT_DIR] [QUERY_PATHS]

Defaults:
  COHORT_DIR  cohort_minsetref_v3
  QUERY_PATHS query_paths.txt

Options:
  --folder DIR              Cohort / graph-build directory (default:
                            cohort_minsetref_v3). Both the Snakemake layout
                            (DIR/Graphs, DIR/inputs/reference_alternatives.fa)
                            and the legacy layout (DIR/localblocks,
                            DIR/main_chroms.fa) are recognised.
  --query FILE              Assembly query-path list (default: query_paths.txt)
  --large-novels FILE       Established large novel loci FASTA; may be empty
                            (default: COHORT_DIR/novel_loci.fa)
  --reference FASTA         Cleaned main reference catalog with source=
                            coordinate headers (default:
                            COHORT_DIR/inputs/reference_alternatives.fa, else
                            COHORT_DIR/main_chroms.fa)
  --reference-haplotype H   Primary reference haplotype name used for the
                            coordinate analysis and as the first small-novel
                            priority tier (default: CHM13_h1)
  --priority LIST           Comma-separated small-novel priority selectors
                            (default: REFERENCE_HAPLOTYPE,HG38_h1,HG002)
  --local-cores N           Total CPU budget (default: 64). Local mapping uses
                            floor(N/4) jobs with 4 cores per job; global
                            alignment and refinement use N cores/jobs.
  --annotations-only        Snakemake mode. Skip the coordinate analysis,
                            unified refinement, and LocalGraphs stages (they
                            are separate workflow rules), run only the
                            local/global/small novel mapping stages, and write
                            immutable sidecar annotations plus
                            _ANNOTATIONS_SUCCESS.json under
                            COHORT_DIR/cohort_annotations/. No graph,
                            alignment, or linear file is rewritten.
  --skip-blastn             Skip the BLASTN residual passes of the local,
                            global, and small-novel mapping stages
  -h, --help                Show this help

Without --annotations-only the seven stages run sequentially and all use
--resume. A failed stage stops the driver; running the same command again
resumes completed work.
EOF
}

total_cores=64
cores_per_local_job=4
folder_option=""
query_option=""
large_novels_option=""
reference_option=""
reference_haplotype="CHM13_h1"
priority_option=""
annotations_only=0
skip_blastn=0
positional=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --folder)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; exit 2; }
            folder_option="$2"
            shift 2
            ;;
        --query)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; exit 2; }
            query_option="$2"
            shift 2
            ;;
        --large-novels)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; exit 2; }
            large_novels_option="$2"
            shift 2
            ;;
        --reference)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; exit 2; }
            reference_option="$2"
            shift 2
            ;;
        --reference-haplotype)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; exit 2; }
            reference_haplotype="$2"
            shift 2
            ;;
        --priority)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; exit 2; }
            priority_option="$2"
            shift 2
            ;;
        --local-cores)
            [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a value" >&2; exit 2; }
            total_cores="$2"
            shift 2
            ;;
        --annotations-only)
            annotations_only=1
            shift
            ;;
        --skip-blastn)
            skip_blastn=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            positional+=("$@")
            break
            ;;
        -*)
            echo "ERROR: unknown option: $1" >&2
            usage
            exit 2
            ;;
        *)
            positional+=("$1")
            shift
            ;;
    esac
done

if [[ ${#positional[@]} -gt 2 ]]; then
    usage
    exit 2
fi
if [[ -n "$folder_option" && ${#positional[@]} -ge 1 ]]; then
    echo "ERROR: specify the cohort with either --folder or a positional argument, not both" >&2
    exit 2
fi
if [[ -n "$query_option" && ${#positional[@]} -ge 2 ]]; then
    echo "ERROR: specify the query list with either --query or a positional argument, not both" >&2
    exit 2
fi
if [[ ! "$total_cores" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --local-cores must be a positive integer; got: $total_cores" >&2
    exit 2
fi
if (( total_cores < cores_per_local_job )); then
    echo "ERROR: --local-cores must be at least $cores_per_local_job" >&2
    exit 2
fi
if [[ -z "$reference_haplotype" ]]; then
    echo "ERROR: --reference-haplotype must not be empty" >&2
    exit 2
fi
local_jobs=$((total_cores / cores_per_local_job))

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
launch_dir="$PWD"
cohort_input="${folder_option:-${positional[0]:-cohort_minsetref_v3}}"
query_input="${query_option:-${positional[1]:-query_paths.txt}}"
python_bin="${MINSETREF_PYTHON:-python3}"

resolve_path() {
    # Absolute paths pass through; relative paths resolve against the launch
    # directory, matching the original driver's behaviour.
    if [[ "$1" == /* ]]; then
        printf '%s\n' "$1"
    else
        printf '%s/%s\n' "$launch_dir" "$1"
    fi
}

cohort_dir="$(resolve_path "$cohort_input")"
query_paths="$(resolve_path "$query_input")"

# Per-partition tree: Snakemake layout first, then the legacy layout.
if [[ -d "$cohort_dir/Graphs" ]]; then
    localblocks="$cohort_dir/Graphs"
else
    localblocks="$cohort_dir/localblocks"
fi

# Cleaned main reference: explicit option, Snakemake catalog, legacy file.
if [[ -n "$reference_option" ]]; then
    main_reference="$(resolve_path "$reference_option")"
elif [[ -s "$cohort_dir/inputs/reference_alternatives.fa" ]]; then
    main_reference="$cohort_dir/inputs/reference_alternatives.fa"
else
    main_reference="$cohort_dir/main_chroms.fa"
fi

coordinate_analysis="$cohort_dir/partition_coordinate_analysis"
unified_refinement="$cohort_dir/unified_novel_refinement"
if (( annotations_only )); then
    # Everything this pass writes lives under one sidecar root so the graph
    # build tree (Graphs/, summary/, checkpoints/) stays untouched.
    annotation_root="$cohort_dir/cohort_annotations"
    partition_novelty="$annotation_root/partition_novelty"
    global_mapping="$annotation_root/global_novel_mapping"
    small_novel_mapping="$annotation_root/small_novel_mapping"
    adjusted_paths="$annotation_root"
    local_graphs=""
else
    annotation_root=""
    partition_novelty="$cohort_dir/partition_novelty"
    global_mapping="$cohort_dir/global_novel_mapping"
    small_novel_mapping="$cohort_dir/small_novel_mapping"
    adjusted_paths="$cohort_dir/adjusted_partition_paths"
    local_graphs="$cohort_dir/LocalGraphs"
fi

if [[ -n "$large_novels_option" ]]; then
    large_novels="$(resolve_path "$large_novels_option")"
else
    large_novels="$cohort_dir/novel_loci.fa"
fi

if [[ -n "$priority_option" ]]; then
    priority="$priority_option"
else
    # Reference haplotype first, then the historical cohort tiers, without
    # repeating the reference when it is one of them.
    priority="$reference_haplotype"
    for tier in HG38_h1 HG002; do
        if [[ "$tier" != "$reference_haplotype" ]]; then
            priority="$priority,$tier"
        fi
    done
fi

if [[ ! -d "$localblocks" ]]; then
    echo "ERROR: partition directory does not exist: $cohort_dir/Graphs (or legacy $cohort_dir/localblocks)" >&2
    exit 1
fi
if [[ ! -s "$query_paths" ]]; then
    echo "ERROR: query-path list is missing or empty: $query_paths" >&2
    exit 1
fi
if [[ ! -s "$main_reference" ]]; then
    echo "ERROR: main reference is missing or empty: $main_reference" >&2
    exit 1
fi
if [[ ! -f "$large_novels" ]]; then
    echo "ERROR: large novel-locus FASTA does not exist: $large_novels" >&2
    exit 1
fi
if (( annotations_only )); then
    # The Snakemake workflow owns these two stages; refuse to silently
    # recompute them with possibly different settings.
    for marker in \
        "$coordinate_analysis/_COORDINATE_ANALYSIS_SUCCESS.json" \
        "$unified_refinement/_UNIFIED_TEMPLATE_REFINEMENT_SUCCESS.json"
    do
        if [[ ! -s "$marker" ]]; then
            echo "ERROR: --annotations-only requires the completed workflow stage marker: $marker" >&2
            exit 1
        fi
    done
fi

skip_blastn_args=()
if (( skip_blastn )); then
    skip_blastn_args=(--skip-blastn)
fi

echo "Cohort directory: $cohort_dir"
echo "Partition directory: $localblocks"
echo "Main reference: $main_reference"
echo "Reference haplotype: $reference_haplotype"
echo "Small-novel priority: $priority"
echo "Large novel loci: $large_novels"
if (( annotations_only )); then
    echo "Mode: annotations only -> $annotation_root"
else
    echo "Mode: full seven-stage run"
fi
if (( skip_blastn )); then
    echo "BLASTN residual passes: skipped"
fi
echo "Total CPU budget: $total_cores"
echo "Analyze/local/global-extraction jobs: $local_jobs"
echo "Coordinate sweep jobs: $local_jobs"
echo "Local alignment: $local_jobs jobs x $cores_per_local_job cores"
echo "Unified nearby-template alignment: $local_jobs jobs x $cores_per_local_job cores"
echo "Global alignment cores: $total_cores"
echo "Refinement jobs: $total_cores"

current_stage="startup"
trap 'rc=$?; echo "ERROR: stage ${current_stage} failed with exit code ${rc}" >&2; exit "$rc"' ERR

run_stage() {
    current_stage="$1"
    shift
    echo
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] START: $current_stage"
    "$@"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE:  $current_stage"
}

if (( ! annotations_only )); then
    run_stage "analyze partition coordinates" \
        "$python_bin" "$script_dir/analyze_partition_coordinates.py" \
            -i "$localblocks" \
            -o "$coordinate_analysis" \
            -j "$local_jobs" \
            --sweep-jobs "$local_jobs" \
            --input-anchor 50 \
            --reference-haplotypes "$reference_haplotype" \
            --ignore-contigs chrM \
            --resume

    run_stage "refine templates from all unsplit novel intervals" \
        "$python_bin" "$script_dir/map_all_novels_to_local_templates.py" \
            -a "$coordinate_analysis" \
            -q "$query_paths" \
            -r "$main_reference" \
            -o "$unified_refinement" \
            -j "$local_jobs" \
            -c "$cores_per_local_job" \
            --minimum-flank 5000 \
            --maximum-template-gap 100 \
            --minimum-unmasked 100 \
            --resume
fi

# ${arr[@]+"${arr[@]}"} expands to nothing for an empty array under `set -u`
# on bash 3.2 (macOS), where a bare "${arr[@]}" is an unbound-variable error.
run_stage "map partition-local novel intervals" \
    "$python_bin" "$script_dir/map_partition_local_novel.py" \
        -a "$coordinate_analysis" \
        -o "$partition_novelty" \
        -q "$query_paths" \
        -j "$local_jobs" \
        -c "$cores_per_local_job" \
        --input-anchor 50 \
        --merge-gap 500 \
        --min-score 100 \
        --min-identity 95 \
        ${skip_blastn_args[@]+"${skip_blastn_args[@]}"} \
        --resume

run_stage "map global novel intervals" \
    "$python_bin" "$script_dir/map_global_novel.py" \
        -a "$coordinate_analysis" \
        -l "$partition_novelty" \
        -r "$main_reference" \
        -q "$query_paths" \
        -o "$global_mapping" \
        -j "$local_jobs" \
        -c "$total_cores" \
        --min-score 100 \
        --min-identity 95 \
        --min-lift-unmasked 100 \
        --minimap-query-batch 50M \
        ${skip_blastn_args[@]+"${skip_blastn_args[@]}"} \
        --resume

run_stage "build and annotate the small-novel minset" \
    "$python_bin" "$script_dir/map_small_novel_minset.py" \
        -a "$coordinate_analysis" \
        -g "$global_mapping" \
        -l "$partition_novelty" \
        -n "$large_novels" \
        -q "$query_paths" \
        -o "$small_novel_mapping" \
        -j "$local_jobs" \
        -c "$total_cores" \
        -p "$priority" \
        --alignment-anchor 50 \
        --min-score 100 \
        --min-identity 95 \
        --minimap-query-batch 50M \
        ${skip_blastn_args[@]+"${skip_blastn_args[@]}"} \
        --resume

annotations_only_args=()
if (( annotations_only )); then
    annotations_only_args=(--annotations-only)
    apply_stage_name="annotate residual novels (immutable sidecars)"
else
    apply_stage_name="apply fixed templates and annotate residual novels"
fi

run_stage "$apply_stage_name" \
    "$python_bin" "$script_dir/apply_unified_novel_refinement.py" \
        -a "$coordinate_analysis" \
        -u "$unified_refinement" \
        -q "$query_paths" \
        -r "$main_reference" \
        -o "$adjusted_paths" \
        --global-queries "$global_mapping/all_global_queries.tsv" \
        --alignment-file "$partition_novelty/all_local_novel_reporting_alignments.tsv" \
        --alignment-file "$global_mapping/all_global_reference_reporting_alignments.tsv" \
        --alignment-file "$small_novel_mapping/all_small_novel_alignments.tsv" \
        --target-fasta "large_novel=$large_novels" \
        --target-fasta "small_novel=$small_novel_mapping/smallnovel.fa" \
        --min-score 100 \
        -j "$total_cores" \
        ${annotations_only_args[@]+"${annotations_only_args[@]}"} \
        --resume

if (( annotations_only )); then
    if [[ ! -s "$annotation_root/_ANNOTATIONS_SUCCESS.json" ]]; then
        echo "ERROR: annotation stage finished without writing $annotation_root/_ANNOTATIONS_SUCCESS.json" >&2
        exit 1
    fi
    current_stage="complete"
    trap - ERR
    echo
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Novel-sequence annotation sidecars completed under $annotation_root."
    exit 0
fi

run_stage "build refined LocalGraphs package" \
    "$python_bin" "$script_dir/build_local_graphs.py" \
        -i "$adjusted_paths" \
        -q "$query_paths" \
        --reference-haplotype "$reference_haplotype" \
        -o "$local_graphs" \
        -j "$total_cores" \
        --resume

current_stage="complete"
trap - ERR
echo
echo "[$(date '+%Y-%m-%d %H:%M:%S')] All seven partition-novelty stages completed."
