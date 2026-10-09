"""
benchmark.py -- Performance benchmark for the authenticated ECIES-like framework.

Paper
-----
A Genus-1 Reformulation of an alpha-beta Modified Edwards-Type Curve
for Secure Telemedicine Data Transmission

Authors: P. Adeola, A. T. Akinwale, B. B. Olofin, A. E. Babalola
         Department of Computing, Anchor University, Lagos, Nigeria

What this script measures
--------------------------
For each of three plaintext sizes (13, 55, and 1000 bytes) and each of
three curve implementations (jq255e, X25519, P-256), it collects 50 paired
timing observations (in microseconds) of:

  - ECDH / shared-secret computation
  - Full encryption  (ephemeral key + ECDH + BLAKE2s KDF + AES-256-GCM)
  - Full decryption  (ECDH + BLAKE2s KDF + AES-256-GCM verify + decrypt)

Key generation (independent of message size) is measured separately and
printed to the console; it is not included in the paired CSV because the
paired comparisons in Table 4 of the manuscript cover ECDH, encryption,
and decryption only.

Both the jq255e sender and recipient keys are generated fresh for each
trial block; the same block index is used for all three curves so that
paired differences eliminate within-block system-load variation.

Output
------
  data/paired_benchmark_raw.csv  -- one row per (curve, operation, trial)

Dependencies
------------
  pip install cryptography
  (jq255e.py must be in the same directory)

Reproducibility note
--------------------
The reported results in the manuscript were obtained on Windows 11,
Intel Core i5, 8 GB RAM, Python 3.12.2.  Execution times will differ
on other platforms because the jq255e construction is implemented in
pure Python affine coordinates, while X25519 and P-256 rely on
optimised C library routines.  See Section 5.4 of the manuscript for
a full discussion of benchmark limitations.
"""

import csv
import hashlib
import os
import secrets
import time
from pathlib import Path

# -- third-party --
from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey
from cryptography.hazmat.primitives.asymmetric.ec import (
    SECP256R1,
    ECDH,
    generate_private_key,
)
from cryptography.hazmat.primitives.ciphers.aead import AESGCM
from cryptography.hazmat.backends import default_backend

# -- local --
from jq255e import H, r, JQ255ePoint, validate_public_key

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

TRIALS      = 50
WARMUP      = 5
MSG_SIZES   = [13, 55, 1000]
OUTPUT_PATH = Path("data") / "paired_benchmark_raw.csv"

# ---------------------------------------------------------------------------
# Utility: timing
# ---------------------------------------------------------------------------

def now_us():
    """Current time in microseconds (high-resolution)."""
    return time.perf_counter_ns() / 1_000.0


# ---------------------------------------------------------------------------
# Utility: BLAKE2s key derivation
# ---------------------------------------------------------------------------

def blake2s_kdf(shared_point_bytes: bytes) -> bytes:
    """Derive a 32-byte AES-256 key from the encoded shared point."""
    return hashlib.blake2s(shared_point_bytes, digest_size=32).digest()


# ---------------------------------------------------------------------------
# Utility: AES-256-GCM
# ---------------------------------------------------------------------------

def aes_gcm_encrypt(key: bytes, plaintext: bytes):
    nonce = os.urandom(12)
    ct    = AESGCM(key).encrypt(nonce, plaintext, None)
    return nonce, ct

def aes_gcm_decrypt(key: bytes, nonce: bytes, ciphertext: bytes):
    return AESGCM(key).decrypt(nonce, ciphertext, None)


# ---------------------------------------------------------------------------
# jq255e ECIES-like operations
# ---------------------------------------------------------------------------

def jq255e_keygen():
    """Generate a jq255e private/public key pair."""
    priv = secrets.randbelow(r - 1) + 1
    pub  = H * priv
    return priv, pub

def jq255e_ecdh(sender_priv: int, recipient_pub: JQ255ePoint) -> bytes:
    """Compute the ECDH shared point and return its encoding."""
    validate_public_key(recipient_pub)
    shared = recipient_pub * sender_priv
    return shared.encode()

def jq255e_encrypt(recipient_pub: JQ255ePoint, plaintext: bytes):
    """
    Full encryption:
      1. Generate ephemeral key pair.
      2. Compute ECDH shared secret.
      3. Derive AES-256 key via BLAKE2s.
      4. Encrypt with AES-256-GCM.
    Returns (ephemeral_pub, nonce, ciphertext).
    """
    ke, ephem_pub = jq255e_keygen()
    shared_bytes  = jq255e_ecdh(ke, recipient_pub)
    session_key   = blake2s_kdf(shared_bytes)
    nonce, ct     = aes_gcm_encrypt(session_key, plaintext)
    return ephem_pub, nonce, ct

def jq255e_decrypt(recipient_priv: int, ephem_pub: JQ255ePoint,
                   nonce: bytes, ciphertext: bytes) -> bytes:
    """
    Full decryption:
      1. Validate ephemeral public key.
      2. Recompute shared secret.
      3. Derive AES-256 key via BLAKE2s.
      4. Authenticate and decrypt.
    """
    shared_bytes = jq255e_ecdh(recipient_priv, ephem_pub)
    session_key  = blake2s_kdf(shared_bytes)
    return aes_gcm_decrypt(session_key, nonce, ciphertext)


# ---------------------------------------------------------------------------
# X25519 ECIES-like operations
# ---------------------------------------------------------------------------

def x25519_keygen():
    priv = X25519PrivateKey.generate()
    pub  = priv.public_key()
    return priv, pub

def x25519_ecdh(sender_priv, recipient_pub) -> bytes:
    return sender_priv.exchange(recipient_pub)

def x25519_encrypt(recipient_pub, plaintext: bytes):
    ke, ephem_pub = x25519_keygen()
    shared_bytes  = x25519_ecdh(ke, recipient_pub)
    session_key   = blake2s_kdf(shared_bytes)
    nonce, ct     = aes_gcm_encrypt(session_key, plaintext)
    return ephem_pub, nonce, ct

def x25519_decrypt(recipient_priv, ephem_pub, nonce: bytes,
                   ciphertext: bytes) -> bytes:
    shared_bytes = x25519_ecdh(recipient_priv, ephem_pub)
    session_key  = blake2s_kdf(shared_bytes)
    return aes_gcm_decrypt(session_key, nonce, ciphertext)


# ---------------------------------------------------------------------------
# P-256 ECIES-like operations
# ---------------------------------------------------------------------------

_BACKEND = default_backend()

def p256_keygen():
    priv = generate_private_key(SECP256R1(), _BACKEND)
    pub  = priv.public_key()
    return priv, pub

def p256_ecdh(sender_priv, recipient_pub) -> bytes:
    return sender_priv.exchange(ECDH(), recipient_pub)

def p256_encrypt(recipient_pub, plaintext: bytes):
    ke, ephem_pub = p256_keygen()
    shared_bytes  = p256_ecdh(ke, recipient_pub)
    session_key   = blake2s_kdf(shared_bytes)
    nonce, ct     = aes_gcm_encrypt(session_key, plaintext)
    return ephem_pub, nonce, ct

def p256_decrypt(recipient_priv, ephem_pub, nonce: bytes,
                 ciphertext: bytes) -> bytes:
    shared_bytes = p256_ecdh(recipient_priv, ephem_pub)
    session_key  = blake2s_kdf(shared_bytes)
    return aes_gcm_decrypt(session_key, nonce, ciphertext)


# ---------------------------------------------------------------------------
# Key-generation benchmark (console output only, not paired)
# ---------------------------------------------------------------------------

def benchmark_keygen(warmup=WARMUP, trials=TRIALS):
    print("\n-- Key Generation (independent of message size) --")
    results = {}
    for name, fn in [("jq255e", jq255e_keygen),
                     ("X25519", x25519_keygen),
                     ("P-256",  p256_keygen)]:
        for _ in range(warmup):
            fn()
        times = []
        for _ in range(trials):
            t0 = now_us()
            fn()
            times.append(now_us() - t0)
        mean = sum(times) / len(times)
        std  = (sum((t - mean) ** 2 for t in times) / (len(times) - 1)) ** 0.5
        results[name] = (mean, std)
        print(f"  {name:8s}: {mean:10.2f} +/- {std:.2f} us  (n={trials})")
    return results


# ---------------------------------------------------------------------------
# Paired benchmark: ECDH, encryption, decryption
# ---------------------------------------------------------------------------

def run_paired_benchmark():
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    rows = []

    for msg_size in MSG_SIZES:
        plaintext = secrets.token_bytes(msg_size)
        print(f"\n=== Message size: {msg_size} bytes ===")

        # -- Warm-up for all three curves --
        print("  Warming up...")
        for _ in range(WARMUP):
            jq_priv_r, jq_pub_r   = jq255e_keygen()
            jq_priv_s, _          = jq255e_keygen()
            jq255e_ecdh(jq_priv_s, jq_pub_r)
            epk, n, ct = jq255e_encrypt(jq_pub_r, plaintext)
            jq255e_decrypt(jq_priv_r, epk, n, ct)

            x25_priv_r, x25_pub_r = x25519_keygen()
            x25_priv_s, _         = x25519_keygen()
            x25519_ecdh(x25_priv_s, x25_pub_r)
            epk, n, ct = x25519_encrypt(x25_pub_r, plaintext)
            x25519_decrypt(x25_priv_r, epk, n, ct)

            p256_priv_r, p256_pub_r = p256_keygen()
            p256_priv_s, _          = p256_keygen()
            p256_ecdh(p256_priv_s, p256_pub_r)
            epk, n, ct = p256_encrypt(p256_pub_r, plaintext)
            p256_decrypt(p256_priv_r, epk, n, ct)

        # -- Paired trials --
        print("  Running paired trials...")
        for trial in range(TRIALS):
            # Fresh keys for this trial block (same block for all curves)
            jq_priv_r, jq_pub_r     = jq255e_keygen()
            jq_priv_s, _            = jq255e_keygen()
            x25_priv_r, x25_pub_r   = x25519_keygen()
            x25_priv_s, _           = x25519_keygen()
            p256_priv_r, p256_pub_r = p256_keygen()
            p256_priv_s, _          = p256_keygen()

            # ---- jq255e ----
            # ECDH
            t0 = now_us()
            shared_bytes_jq = jq255e_ecdh(jq_priv_s, jq_pub_r)
            rows.append((msg_size, "ecdh", "jq255e", trial, now_us() - t0))

            # Encryption
            t0 = now_us()
            epk_jq, nonce_jq, ct_jq = jq255e_encrypt(jq_pub_r, plaintext)
            rows.append((msg_size, "encryption", "jq255e", trial, now_us() - t0))

            # Decryption
            t0 = now_us()
            jq255e_decrypt(jq_priv_r, epk_jq, nonce_jq, ct_jq)
            rows.append((msg_size, "decryption", "jq255e", trial, now_us() - t0))

            # ---- X25519 ----
            t0 = now_us()
            x25519_ecdh(x25_priv_s, x25_pub_r)
            rows.append((msg_size, "ecdh", "X25519", trial, now_us() - t0))

            t0 = now_us()
            epk_x25, nonce_x25, ct_x25 = x25519_encrypt(x25_pub_r, plaintext)
            rows.append((msg_size, "encryption", "X25519", trial, now_us() - t0))

            t0 = now_us()
            x25519_decrypt(x25_priv_r, epk_x25, nonce_x25, ct_x25)
            rows.append((msg_size, "decryption", "X25519", trial, now_us() - t0))

            # ---- P-256 ----
            t0 = now_us()
            p256_ecdh(p256_priv_s, p256_pub_r)
            rows.append((msg_size, "ecdh", "P-256", trial, now_us() - t0))

            t0 = now_us()
            epk_p256, nonce_p256, ct_p256 = p256_encrypt(p256_pub_r, plaintext)
            rows.append((msg_size, "encryption", "P-256", trial, now_us() - t0))

            t0 = now_us()
            p256_decrypt(p256_priv_r, epk_p256, nonce_p256, ct_p256)
            rows.append((msg_size, "decryption", "P-256", trial, now_us() - t0))

            if (trial + 1) % 10 == 0:
                print(f"    Trial {trial + 1}/{TRIALS} complete.")

    # -- Write CSV --
    with open(OUTPUT_PATH, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["message_size_bytes", "operation", "curve",
                         "trial", "time_us"])
        writer.writerows(rows)

    print(f"\nResults saved to {OUTPUT_PATH}")
    return rows


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("jq255e ECIES-like Framework -- Performance Benchmark")
    print("=" * 60)
    print(f"Trials per condition : {TRIALS}")
    print(f"Warm-up iterations   : {WARMUP}")
    print(f"Message sizes (bytes): {MSG_SIZES}")
    print(f"Curves               : jq255e, X25519, P-256")

    benchmark_keygen()
    run_paired_benchmark()

    print("\nDone.")
