#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLIENT="$SCRIPT_DIR/threshold-client"
SERVER="$SCRIPT_DIR/threshold-server"

TAG_HEX="${TAG_HEX:-504f5052462d7461672d303030303100}"
X_HEX="${X_HEX:-5ac31980277104bed261}"
NSERVERS="${NSERVERS:-3}"
THRESHOLD="${THRESHOLD:-2}"

default_q_spec() {
  local threshold="$1"
  local i
  local out=""

  for ((i = 0; i < threshold; i++)); do
    if [[ -n "$out" ]]; then
      out+=","
    fi
    out+="$i"
  done
  printf '%s' "$out"
}

Q_SPEC="${Q_SPEC:-$(default_q_spec "$THRESHOLD")}"

get_field() {
  local key="$1"
  awk -F= -v k="$key" '$1 == k { print $2 }'
}

add_fp() {
  awk -v a="$1" -v b="$2" 'BEGIN { printf "%.6f", a + b }'
}

mul_fp() {
  awk -v a="$1" -v b="$2" 'BEGIN { printf "%.6f", a * b }'
}

sub_fp() {
  awk -v a="$1" -v b="$2" 'BEGIN { printf "%.6f", a - b }'
}

bytes_to_kb() {
  awk -v b="$1" 'BEGIN { printf "%.2f", b / 1024 }'
}

require_ok() {
  local output="$1"
  local label="$2"
  local ok
  ok="$(printf '%s\n' "$output" | get_field ok)"
  if [[ "$ok" != "1" ]]; then
    printf 'step failed: %s\n' "$label" >&2
    printf '%s\n' "$output" >&2
    exit 1
  fi
}

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

pp_file="$tmpdir/pp.bin"
master_file="$tmpdir/master.bin"
pk_dir="$tmpdir/pk"
sk_dir="$tmpdir/sk"
preq_file="$tmpdir/preq.bin"
pre_state_file="$tmpdir/pre_state.bin"
prep_dir="$tmpdir/prep"
server_pre_dir="$tmpdir/server_pre"
req_file="$tmpdir/req.bin"
st_file="$tmpdir/st.bin"
rep_dir="$tmpdir/rep"
SLOTS="${SLOTS:-10}"

setup_out="$("$CLIENT" setup "$pp_file")"
require_ok "$setup_out" "setup"

keygen_out="$("$SERVER" keygen "$pp_file" "$NSERVERS" "$THRESHOLD" "$pk_dir" "$sk_dir" "$master_file")"
require_ok "$keygen_out" "keygen"

preproc_client_out="$("$CLIENT" preproc-client "$pp_file" "$SLOTS" "$preq_file" "$pre_state_file")"
require_ok "$preproc_client_out" "preproc-client"

preproc_server_out="$("$SERVER" preproc-server "$pp_file" "$sk_dir" "$NSERVERS" "$THRESHOLD" "$preq_file" "$prep_dir" "$server_pre_dir")"
require_ok "$preproc_server_out" "preproc-server"

request_out="$("$CLIENT" request "$pp_file" "$pre_state_file" "$NSERVERS" "$THRESHOLD" "$Q_SPEC" "$TAG_HEX" "$X_HEX" "$req_file" "$st_file")"
require_ok "$request_out" "request"

if [[ -n "${LEOPARD_REQUEST_RAW_LOG:-}" ]]; then
  printf '%s\n' "$request_out" > "$LEOPARD_REQUEST_RAW_LOG"
fi
# Each online server checks the request independently; do not extrapolate a
# single verifier measurement. Both verification and proof jobs are sequential.
IFS=',' read -r -a online_servers <<< "$Q_SPEC"
for server in "${online_servers[@]}"; do
  one_verify="$("$CLIENT" verify-request "$pp_file" "$TAG_HEX" "$req_file")"
  require_ok "$one_verify" "verify-request server $server"
  printf '%s\n' "$one_verify" >> "$tmpdir/request_verifiers.log"
done
verify_req_out="$(awk -F= '/^[a-zA-Z0-9_]+=[0-9.]+$/ {sum[$1]+=$2;n[$1]++}
  END {for(k in sum) printf "%s=%.9f\n",k,sum[k]/n[k]}' "$tmpdir/request_verifiers.log")"
blind_verify_total_sec="$(awk -F= '$1=="time_sec" {s+=$2} END {printf "%.9f",s}' "$tmpdir/request_verifiers.log")"

blindeval_core_out="$("$SERVER" blindeval "$pp_file" "$pk_dir" "$sk_dir" "$server_pre_dir" "$TAG_HEX" "$req_file" "$rep_dir")"
require_ok "$blindeval_core_out" "blindeval"

verify_rep_out="$("$SERVER" verify-replies "$pp_file" "$pk_dir" "$prep_dir" "$TAG_HEX" "$req_file" "$rep_dir")"
require_ok "$verify_rep_out" "verify-replies"

finalize_raw_out="$("$CLIENT" finalize "$st_file" "$prep_dir" "$rep_dir")"
require_ok "$finalize_raw_out" "finalize"

eval_out="$("$CLIENT" eval "$pp_file" "$TAG_HEX" "$X_HEX" "$master_file")"
require_ok "$eval_out" "eval"

setup_sec="$(printf '%s\n' "$setup_out" | get_field time_sec)"
keygen_sec="$(printf '%s\n' "$keygen_out" | get_field time_sec)"
keygen_party_total_sec="$(printf '%s\n' "$keygen_out" | get_field keygen_party_total_sec)"
keygen_party_avg_sec="$(printf '%s\n' "$keygen_out" | get_field keygen_party_avg_sec)"
preproc_client_sec="$(printf '%s\n' "$preproc_client_out" | get_field time_sec)"
preproc_server_sec="$(printf '%s\n' "$preproc_server_out" | get_field time_sec)"
preproc_party_total_sec="$(printf '%s\n' "$preproc_server_out" | get_field preproc_party_total_sec)"
preproc_party_avg_sec="$(printf '%s\n' "$preproc_server_out" | get_field preproc_party_avg_sec)"
request_sec="$(printf '%s\n' "$request_out" | get_field time_sec)"
request_slot_j="$(printf '%s\n' "$request_out" | get_field slot_j)"
request_q_mask="$(printf '%s\n' "$request_out" | get_field q_mask)"
req_inner_instance_sec="$(printf '%s\n' "$request_out" | get_field req_inner_instance_sec)"
req_inner_prepare_sec="$(printf '%s\n' "$request_out" | get_field req_inner_prepare_sec)"
req_inner_prove_sec="$(printf '%s\n' "$request_out" | get_field req_inner_prove_sec)"
req_inner_selfcheck_sec="$(printf '%s\n' "$request_out" | get_field req_inner_selfcheck_sec)"
req_inner_rejection_rounds="$(printf '%s\n' "$request_out" | get_field req_inner_rejection_rounds)"
req_inner_batch_rounds="$(printf '%s\n' "$request_out" | get_field req_inner_batch_rounds)"
req_a1open_cache="$(printf '%s\n' "$request_out" | get_field req_a1open_cache)"
req_outer_witness_sec="$(printf '%s\n' "$request_out" | get_field req_outer_witness_sec)"
req_outer_statement_sec="$(printf '%s\n' "$request_out" | get_field req_outer_statement_sec)"
req_outer_stmt_fold_sec="$(printf '%s\n' "$request_out" | get_field req_outer_stmt_fold_sec)"
req_outer_stmt_quad_sec="$(printf '%s\n' "$request_out" | get_field req_outer_stmt_quad_sec)"
req_outer_stmt_cnst1_sec="$(printf '%s\n' "$request_out" | get_field req_outer_stmt_cnst1_sec)"
req_outer_stmt_cnst2_sec="$(printf '%s\n' "$request_out" | get_field req_outer_stmt_cnst2_sec)"
req_outer_stmt_cnst3_sec="$(printf '%s\n' "$request_out" | get_field req_outer_stmt_cnst3_sec)"
req_outer_stmt_hash_sec="$(printf '%s\n' "$request_out" | get_field req_outer_stmt_hash_sec)"
req_outer_prove_sec="$(printf '%s\n' "$request_out" | get_field req_outer_prove_sec)"
req_serialize_sec="$(printf '%s\n' "$request_out" | get_field req_serialize_sec)"
req_export_sec="$(printf '%s\n' "$request_out" | get_field req_export_sec)"
blind_verify_sec="$(printf '%s\n' "$verify_req_out" | get_field time_sec)"
verify_public_ctx_sec="$(printf '%s\n' "$verify_req_out" | get_field verify_public_ctx_sec)"
verify_a1open_cache="$(printf '%s\n' "$verify_req_out" | get_field verify_a1open_cache)"
verify_fs_seed_sec="$(printf '%s\n' "$verify_req_out" | get_field verify_fs_seed_sec)"
verify_deserialize_sec="$(printf '%s\n' "$verify_req_out" | get_field verify_deserialize_sec)"
verify_outer_statement_sec="$(printf '%s\n' "$verify_req_out" | get_field verify_outer_statement_sec)"
verify_outer_stmt_fold_sec="$(printf '%s\n' "$verify_req_out" | get_field verify_outer_stmt_fold_sec)"
verify_outer_stmt_quad_sec="$(printf '%s\n' "$verify_req_out" | get_field verify_outer_stmt_quad_sec)"
verify_outer_stmt_cnst1_sec="$(printf '%s\n' "$verify_req_out" | get_field verify_outer_stmt_cnst1_sec)"
verify_outer_stmt_cnst2_sec="$(printf '%s\n' "$verify_req_out" | get_field verify_outer_stmt_cnst2_sec)"
verify_outer_stmt_cnst3_sec="$(printf '%s\n' "$verify_req_out" | get_field verify_outer_stmt_cnst3_sec)"
verify_outer_stmt_hash_sec="$(printf '%s\n' "$verify_req_out" | get_field verify_outer_stmt_hash_sec)"
verify_outer_verify_sec="$(printf '%s\n' "$verify_req_out" | get_field verify_outer_verify_sec)"
blind_core_sec="$(printf '%s\n' "$blindeval_core_out" | get_field time_sec)"
blind_q_mask="$(printf '%s\n' "$blindeval_core_out" | get_field q_mask)"
blindeval_bundle_count="$(printf '%s\n' "$blindeval_core_out" | get_field bundle_count)"
blindeval_component_count="$(printf '%s\n' "$blindeval_core_out" | get_field component_count)"
blindeval_build_instance_total_sec="$(printf '%s\n' "$blindeval_core_out" | get_field blindeval_build_instance_total_sec)"
blindeval_prover_total_sec="$(printf '%s\n' "$blindeval_core_out" | get_field blindeval_prover_total_sec)"
blindeval_export_total_sec="$(printf '%s\n' "$blindeval_core_out" | get_field blindeval_export_total_sec)"
blindeval_build_instance_avg_sec="$(printf '%s\n' "$blindeval_core_out" | get_field blindeval_build_instance_avg_sec)"
blindeval_prover_avg_sec="$(printf '%s\n' "$blindeval_core_out" | get_field blindeval_prover_avg_sec)"
blindeval_export_avg_sec="$(printf '%s\n' "$blindeval_core_out" | get_field blindeval_export_avg_sec)"
blind_total_sec="$(add_fp "$blind_verify_total_sec" "$blind_core_sec")"
finalize_verify_sec="$(printf '%s\n' "$verify_rep_out" | get_field time_sec)"
finalize_core_sec="$(printf '%s\n' "$finalize_raw_out" | get_field time_sec)"
finalize_total_sec="$(add_fp "$finalize_verify_sec" "$finalize_core_sec")"
finalize_parallel_verify="$(printf '%s\n' "$verify_rep_out" | get_field parallel_verify_replies)"
finalize_verify_import_total_sec="$(printf '%s\n' "$verify_rep_out" | get_field verify_reply_import_total_sec)"
finalize_verify_build_total_sec="$(printf '%s\n' "$verify_rep_out" | get_field verify_reply_build_instance_total_sec)"
finalize_verify_core_total_sec="$(printf '%s\n' "$verify_rep_out" | get_field verify_reply_verifier_total_sec)"
finalize_verify_profile_total_sec="$(printf '%s\n' "$verify_rep_out" | get_field verify_reply_profile_total_sec)"
finalize_verify_import_avg_sec="$(printf '%s\n' "$verify_rep_out" | get_field verify_reply_import_avg_sec)"
finalize_verify_build_avg_sec="$(printf '%s\n' "$verify_rep_out" | get_field verify_reply_build_instance_avg_sec)"
finalize_verify_core_avg_sec="$(printf '%s\n' "$verify_rep_out" | get_field verify_reply_verifier_avg_sec)"
finalize_verify_profile_avg_sec="$(printf '%s\n' "$verify_rep_out" | get_field verify_reply_profile_avg_sec)"
eval_sec="$(printf '%s\n' "$eval_out" | get_field time_sec)"

pp_bytes="$(printf '%s\n' "$setup_out" | get_field pp_bytes)"
pk_bytes="$(printf '%s\n' "$keygen_out" | get_field pk_total_bytes)"
pk_avg_bytes="$(printf '%s\n' "$keygen_out" | get_field pk_avg_bytes)"
preq_bytes="$(printf '%s\n' "$preproc_client_out" | get_field preq_bytes)"
prep_bytes="$(printf '%s\n' "$preproc_server_out" | get_field prep_total_bytes)"
prep_avg_bytes="$(printf '%s\n' "$preproc_server_out" | get_field prep_avg_bytes)"
req_bytes="$(printf '%s\n' "$request_out" | get_field req_bytes)"
req_public_bytes="$(printf '%s\n' "$request_out" | get_field req_public_bytes)"
req_public_cx_bytes="$(printf '%s\n' "$request_out" | get_field req_public_cx_bytes)"
req_public_other_bytes="$(printf '%s\n' "$request_out" | get_field req_public_other_bytes)"
req_proof_bytes="$(printf '%s\n' "$request_out" | get_field req_proof_bytes)"
req_kb="$(bytes_to_kb "$req_bytes")"
req_public_kb="$(bytes_to_kb "$req_public_bytes")"
req_public_cx_kb="$(bytes_to_kb "$req_public_cx_bytes")"
req_public_other_kb="$(bytes_to_kb "$req_public_other_bytes")"
req_proof_kb="$(bytes_to_kb "$req_proof_bytes")"
req_protocol_input_sec="$(printf '%s\n' "$request_out" | get_field req_protocol_input_sec)"
req_protocol_transcript_sec="$(printf '%s\n' "$request_out" | get_field req_protocol_transcript_sec)"
req_protocol_selector_material_sec="$(printf '%s\n' "$request_out" | get_field req_protocol_selector_material_sec)"
req_protocol_slot_load_sec="$(printf '%s\n' "$request_out" | get_field req_protocol_slot_load_sec)"
req_protocol_hash_eval_sec="$(printf '%s\n' "$request_out" | get_field req_protocol_hash_eval_sec)"
req_protocol_cx_sec="$(printf '%s\n' "$request_out" | get_field req_protocol_cx_sec)"
req_protocol_state_sec="$(printf '%s\n' "$request_out" | get_field req_protocol_state_sec)"
req_protocol_public_export_sec="$(printf '%s\n' "$request_out" | get_field req_protocol_public_export_sec)"
req_fast_public_total_sec="$(printf '%s\n' "$request_out" | get_field req_fast_public_total_sec)"
req_fast_public_selector_material_sec="$(printf '%s\n' "$request_out" | get_field req_fast_public_selector_material_sec)"
req_fast_public_slot_load_sec="$(printf '%s\n' "$request_out" | get_field req_fast_public_slot_load_sec)"
req_fast_public_chain_sec="$(printf '%s\n' "$request_out" | get_field req_fast_public_chain_sec)"
req_fast_public_cx_sec="$(printf '%s\n' "$request_out" | get_field req_fast_public_cx_sec)"
req_nizk_prepare_sec="$(printf '%s\n' "$request_out" | get_field req_nizk_prepare_sec)"
req_nizk_inner_setup_sec="$(printf '%s\n' "$request_out" | get_field req_nizk_inner_setup_sec)"
rep_bytes="$(printf '%s\n' "$blindeval_core_out" | get_field rep_total_bytes)"
rep_avg_bytes="$(printf '%s\n' "$blindeval_core_out" | get_field rep_avg_bytes)"
rep_public_bytes="$(printf '%s\n' "$blindeval_core_out" | get_field rep_public_total_bytes)"
rep_proof_bytes="$(printf '%s\n' "$blindeval_core_out" | get_field rep_proof_total_bytes)"
rep_public_avg_bytes="$(printf '%s\n' "$blindeval_core_out" | get_field rep_public_avg_bytes)"
rep_proof_avg_bytes="$(printf '%s\n' "$blindeval_core_out" | get_field rep_proof_avg_bytes)"
blindeval_protocol_eval_total_sec="$(printf '%s\n' "$blindeval_core_out" | get_field blindeval_protocol_eval_total_sec)"
blindeval_protocol_eval_avg_sec="$(printf '%s\n' "$blindeval_core_out" | get_field blindeval_protocol_eval_avg_sec)"
blindeval_nizk_build_total_sec="$(printf '%s\n' "$blindeval_core_out" | get_field blindeval_nizk_build_total_sec)"
blindeval_nizk_build_avg_sec="$(printf '%s\n' "$blindeval_core_out" | get_field blindeval_nizk_build_avg_sec)"

request_protocol_total_sec="$(add_fp "$req_protocol_input_sec" "$req_protocol_transcript_sec")"
request_protocol_total_sec="$(add_fp "$request_protocol_total_sec" "$req_protocol_state_sec")"
request_protocol_total_sec="$(add_fp "$request_protocol_total_sec" "$req_protocol_public_export_sec")"
request_nizk_total_sec="$(add_fp "$req_nizk_prepare_sec" "$req_nizk_inner_setup_sec")"
request_nizk_total_sec="$(add_fp "$request_nizk_total_sec" "$req_inner_prove_sec")"
request_nizk_total_sec="$(add_fp "$request_nizk_total_sec" "$req_outer_witness_sec")"
request_nizk_total_sec="$(add_fp "$request_nizk_total_sec" "$req_outer_statement_sec")"
request_nizk_total_sec="$(add_fp "$request_nizk_total_sec" "$req_outer_prove_sec")"
request_nizk_total_sec="$(add_fp "$request_nizk_total_sec" "$req_serialize_sec")"
blind_protocol_total_sec="$(add_fp "$blindeval_protocol_eval_total_sec" "$blindeval_export_total_sec")"
blind_server_nizk_wall_sec="$(sub_fp "$blind_core_sec" "$blind_protocol_total_sec")"
blind_nizk_wall_sec="$(add_fp "$blind_verify_total_sec" "$blind_server_nizk_wall_sec")"
blind_nizk_work_total_sec="$(add_fp "$blind_verify_total_sec" "$blindeval_nizk_build_total_sec")"
blind_nizk_work_total_sec="$(add_fp "$blind_nizk_work_total_sec" "$blindeval_prover_total_sec")"

y_finalize="$(printf '%s\n' "$finalize_raw_out" | get_field y_hex)"
y_eval="$(printf '%s\n' "$eval_out" | get_field y_hex)"

correctness="FAIL"
if [[ "$y_finalize" == "$y_eval" ]]; then
  correctness="PASS"
fi

printf 'T-LeOPaRd benchmark\n\n'
printf 'Status:\n'
printf '  correctness       = %s\n' "$correctness"
printf '  threshold mode    = %s-of-%s replicated sharing\n' "$THRESHOLD" "$NSERVERS"
printf '  key commitment    = %s\n' "$(printf '%s\n' "$keygen_out" | get_field commitment_mode)"
printf '  quorum Q          = %s (mask %s)\n' "$Q_SPEC" "${request_q_mask:-$blind_q_mask}"
printf '  preprocessing     = %s client slots, per-server offline vk tables\n' "$SLOTS"
printf '  request proof     = third_party/nizk1 RBIT preprocessed client proof\n'
printf '  server proof      = one third_party/nizk2 same-modulus proof per owner bundle\n'
printf '  reply bundles     = %s bundles covering %s replicated components\n' "$blindeval_bundle_count" "$blindeval_component_count"
printf '  split execution   = client/server tools + driver script\n\n'

printf 'Timing:\n'
printf '  F.Setup           = %s s\n' "$setup_sec"
printf '  F.KeyGen          = %s s\n' "$keygen_sec"
printf '  F.PreProcClient   = %s s\n' "$preproc_client_sec"
printf '  F.PreProcServer   = %s s\n' "$preproc_server_sec"
printf '  F.Request         = %s s\n' "$request_sec"
printf '  F.BlindEval       = %s s\n' "$blind_total_sec"
printf '  F.Finalize        = %s s\n' "$finalize_total_sec"
printf '  F.Eval            = %s s\n\n' "$eval_sec"

printf 'Finalize Breakdown:\n'
printf '  verify replies    = %s s wall\n' "$finalize_verify_sec"
printf '    import u/vk     = %s s total / %s s avg\n' "$finalize_verify_import_total_sec" "$finalize_verify_import_avg_sec"
printf '    build statement = %s s total / %s s avg\n' "$finalize_verify_build_total_sec" "$finalize_verify_build_avg_sec"
printf '    NIZK2 verifier  = %s s total / %s s avg\n' "$finalize_verify_core_total_sec" "$finalize_verify_core_avg_sec"
printf '    profile total   = %s s total / %s s avg\n' "$finalize_verify_profile_total_sec" "$finalize_verify_profile_avg_sec"
printf '  raw finalize      = %s s\n' "$finalize_core_sec"
printf '  parallel verify   = %s\n\n' "${finalize_parallel_verify:-0}"

printf 'KeyGen Breakdown:\n'
printf '  party total       = %s s\n' "$keygen_party_total_sec"
printf '  party avg         = %s s\n\n' "$keygen_party_avg_sec"

printf 'PreProcServer Breakdown:\n'
printf '  party total       = %s s\n' "$preproc_party_total_sec"
printf '  party avg         = %s s\n\n' "$preproc_party_avg_sec"

printf 'Request Breakdown:\n'
printf '  slot J            = %s\n' "$request_slot_j"
printf '  fast public path* = %s s\n' "$req_fast_public_total_sec"
printf '    selector mats   = %s s\n' "$req_fast_public_selector_material_sec"
printf '    load slot       = %s s\n' "$req_fast_public_slot_load_sec"
printf '    x -> b0 chain   = %s s\n' "$req_fast_public_chain_sec"
printf '    build Cx        = %s s\n' "$req_fast_public_cx_sec"
printf '  protocol total    = %s s\n' "$request_protocol_total_sec"
printf '    input hash/bits = %s s\n' "$req_protocol_input_sec"
printf '    selector mats   = %s s\n' "$req_protocol_selector_material_sec"
printf '    load slot       = %s s\n' "$req_protocol_slot_load_sec"
printf '    build B=H(x,t)  = %s s\n' "$req_protocol_hash_eval_sec"
printf '    build Cx        = %s s\n' "$req_protocol_cx_sec"
printf '    state update    = %s s\n' "$req_protocol_state_sec"
printf '    public export   = %s s\n' "$req_protocol_public_export_sec"
printf '  NIZK total        = %s s\n' "$request_nizk_total_sec"
printf '    prepare stmt    = %s s\n' "$req_nizk_prepare_sec"
printf '    inner setup     = %s s\n' "$req_nizk_inner_setup_sec"
printf '  inner instance    = %s s\n' "$req_inner_instance_sec"
printf '  inner prepare     = %s s\n' "$req_inner_prepare_sec"
printf '  inner prove       = %s s\n' "$req_inner_prove_sec"
printf '  inner self-check  = %s s\n' "$req_inner_selfcheck_sec"
printf '  inner reject rnds = %s\n' "$req_inner_rejection_rounds"
printf '  inner batch rnds  = %s\n' "$req_inner_batch_rounds"
printf '  A1open cache      = %s\n' "$req_a1open_cache"
printf '  outer witness     = %s s\n' "$req_outer_witness_sec"
printf '  outer statement   = %s s\n' "$req_outer_statement_sec"
printf '    fold/open       = %s s\n' "$req_outer_stmt_fold_sec"
printf '    quad coeff      = %s s\n' "$req_outer_stmt_quad_sec"
printf '    constraint 1    = %s s\n' "$req_outer_stmt_cnst1_sec"
printf '    constraint 2    = %s s\n' "$req_outer_stmt_cnst2_sec"
printf '    constraint 3    = %s s\n' "$req_outer_stmt_cnst3_sec"
printf '    final hash      = %s s\n' "$req_outer_stmt_hash_sec"
printf '  outer prove       = %s s\n' "$req_outer_prove_sec"
printf '  serialize proof   = %s s\n' "$req_serialize_sec"
printf '  export transcript = %s s\n\n' "$req_export_sec"

printf 'BlindEval Breakdown:\n'
printf '  protocol model    = %s independent request verifications + %s reply bundles\n' "$THRESHOLD" "$blindeval_bundle_count"
printf '  protocol total    = %s s\n' "$blind_protocol_total_sec"
printf '    compute u       = %s s total / %s s avg\n' "$blindeval_protocol_eval_total_sec" "$blindeval_protocol_eval_avg_sec"
printf '    export reply    = %s s total / %s s avg\n' "$blindeval_export_total_sec" "$blindeval_export_avg_sec"
printf '  NIZK wall total   = %s s\n' "$blind_nizk_wall_sec"
printf '  NIZK work total   = %s s\n' "$blind_nizk_work_total_sec"
printf '  verify-request    = %s s per server / %s s total\n' "$blind_verify_sec" "$blind_verify_total_sec"
printf '  A1open cache      = %s\n' "$verify_a1open_cache"
printf '  public ctx        = %s s\n' "$verify_public_ctx_sec"
printf '  FS + ar-seed      = %s s\n' "$verify_fs_seed_sec"
printf '  deserialize req   = %s s\n' "$verify_deserialize_sec"
printf '  outer statement   = %s s\n' "$verify_outer_statement_sec"
printf '    fold/open       = %s s\n' "$verify_outer_stmt_fold_sec"
printf '    quad coeff      = %s s\n' "$verify_outer_stmt_quad_sec"
printf '    constraint 1    = %s s\n' "$verify_outer_stmt_cnst1_sec"
printf '    constraint 2    = %s s\n' "$verify_outer_stmt_cnst2_sec"
printf '    constraint 3    = %s s\n' "$verify_outer_stmt_cnst3_sec"
printf '    final hash      = %s s\n' "$verify_outer_stmt_hash_sec"
printf '  outer verify      = %s s\n' "$verify_outer_verify_sec"
printf '  NIZK2 core total  = %s s\n' "$blind_core_sec"
printf '  NIZK2 build(stmt) = %s s total / %s s avg\n' "$blindeval_nizk_build_total_sec" "$blindeval_nizk_build_avg_sec"
printf '  build instance    = %s s total / %s s avg\n' "$blindeval_build_instance_total_sec" "$blindeval_build_instance_avg_sec"
printf '  NIZK2 prover      = %s s total / %s s avg\n' "$blindeval_prover_total_sec" "$blindeval_prover_avg_sec"
printf '  export reply      = %s s total / %s s avg\n\n' "$blindeval_export_total_sec" "$blindeval_export_avg_sec"

printf 'Communication:\n'
printf '  pp bytes          = %s\n' "$pp_bytes"
printf '  pk bytes          = %s\n' "$pk_bytes"
printf '  pk avg bytes      = %s\n' "$pk_avg_bytes"
printf '  preq bytes        = %s\n' "$preq_bytes"
printf '  prep bytes        = %s\n' "$prep_bytes"
printf '  prep avg bytes    = %s\n' "$prep_avg_bytes"
printf '  req bytes         = %s B = %s KB\n' "$req_bytes" "$req_kb"
printf '  req public bytes  = %s B = %s KB\n' "$req_public_bytes" "$req_public_kb"
printf '    Cx bytes        = %s B = %s KB\n' "$req_public_cx_bytes" "$req_public_cx_kb"
printf '    other public    = %s B = %s KB\n' "$req_public_other_bytes" "$req_public_other_kb"
printf '  req proof bytes   = %s B = %s KB\n' "$req_proof_bytes" "$req_proof_kb"
printf '  rep bytes         = %s\n' "$rep_bytes"
printf '  rep public bytes  = %s\n' "$rep_public_bytes"
printf '  rep proof bytes   = %s\n' "$rep_proof_bytes"
printf '  rep avg bytes     = %s\n\n' "$rep_avg_bytes"

printf 'Outputs:\n'
printf '  Finalize y        = %s\n' "$y_finalize"
printf '  Eval y            = %s\n' "$y_eval"
[[ "$correctness" == PASS ]]
