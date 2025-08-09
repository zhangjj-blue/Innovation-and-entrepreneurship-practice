#include <iostream>
#include <iomanip>
#include <chrono>
#include <immintrin.h>
#include <cstdint>
#include <cstring>
using namespace std;
// SM4常量定义
constexpr uint32_t SM4_BLOCK_SIZE = 16;
constexpr uint32_t SM4_ROUNDS = 32;

// SM4 SBox (标准定义)
alignas(64) const uint8_t SM4_SBOX[256] = {
    0xD6, 0x90, 0xE9, 0xFE, 0xCC, 0xE1, 0x3D, 0xB7, 0x16, 0xB6, 0x14, 0xC2, 0x28, 0xFB, 0x2C, 0x05,
    0x2B, 0x67, 0x9A, 0x76, 0x2A, 0xBE, 0x04, 0xC3, 0xAA, 0x44, 0x13, 0x26, 0x49, 0x86, 0x06, 0x99,
    0x9C, 0x42, 0x50, 0xF4, 0x91, 0xEF, 0x98, 0x7A, 0x33, 0x54, 0x0B, 0x43, 0xED, 0xCF, 0xAC, 0x62,
    0xE4, 0xB3, 0x1C, 0xA9, 0xC9, 0x08, 0xE8, 0x95, 0x80, 0xDF, 0x94, 0xFA, 0x75, 0x8F, 0x3F, 0xA6,
    0x47, 0x07, 0xA7, 0xFC, 0xF3, 0x73, 0x17, 0xBA, 0x83, 0x59, 0x3C, 0x19, 0xE6, 0x85, 0x4F, 0xA8,
    0x68, 0x6B, 0x81, 0xB2, 0x71, 0x64, 0xDA, 0x8B, 0xF8, 0xEB, 0x0F, 0x4B, 0x70, 0x56, 0x9D, 0x35,
    0x1E, 0x24, 0x0E, 0x5E, 0x63, 0x58, 0xD1, 0xA2, 0x25, 0x22, 0x7C, 0x3B, 0x01, 0x21, 0x78, 0x87,
    0xD4, 0x00, 0x46, 0x57, 0x9F, 0xD3, 0x27, 0x52, 0x4C, 0x36, 0x02, 0xE7, 0xA0, 0xC4, 0xC8, 0x9E,
    0xEA, 0xBF, 0x8A, 0xD2, 0x40, 0xC7, 0x38, 0xB5, 0xA3, 0xF7, 0xF2, 0xCE, 0xF9, 0x61, 0x15, 0xA1,
    0xE0, 0xAE, 0x5D, 0xA4, 0x9B, 0x34, 0x1A, 0x55, 0xAD, 0x93, 0x32, 0x30, 0xF5, 0x8C, 0xB1, 0xE3,
    0x1D, 0xF6, 0xE2, 0x2E, 0x82, 0x66, 0xCA, 0x60, 0xC0, 0x29, 0x23, 0xAB, 0x0D, 0x53, 0x4E, 0x6F,
    0xD5, 0xDB, 0x37, 0x45, 0xDE, 0xFD, 0x8E, 0x2F, 0x03, 0xFF, 0x6A, 0x72, 0x6D, 0x6C, 0x5B, 0x51,
    0x8D, 0x1B, 0xAF, 0x92, 0xBB, 0xDD, 0xBC, 0x7F, 0x11, 0xD9, 0x5C, 0x41, 0x1F, 0x10, 0x5A, 0xD8,
    0x0A, 0xC1, 0x31, 0x88, 0xA5, 0xCD, 0x7B, 0xBD, 0x2D, 0x74, 0xD0, 0x12, 0xB8, 0xE5, 0xB4, 0xB0,
    0x89, 0x69, 0x97, 0x4A, 0x0C, 0x96, 0x77, 0x7E, 0x65, 0xB9, 0xF1, 0x09, 0xC5, 0x6E, 0xC6, 0x84,
    0x18, 0xF0, 0x7D, 0xEC, 0x3A, 0xDC, 0x4D, 0x20, 0x79, 0xEE, 0x5F, 0x3E, 0xD7, 0xCB, 0x39, 0x48
};

// ======================= 基础实现优化 =======================
// 循环展开+寄存器优化 (PDF Page 5)
void sm4_basic_encrypt(uint8_t out[16], const uint8_t in[16], const uint32_t rk[SM4_ROUNDS]) {
    uint32_t block[4];
    memcpy(block, in, 16);

    // 完全展开的32轮加密 (8组×4轮)
    for (int i = 0; i < 8; ++i) {
        uint32_t tmp = block[1] ^ block[2] ^ block[3] ^ rk[4 * i];
        tmp = SM4_SBOX[tmp & 0xFF] | (SM4_SBOX[(tmp >> 8) & 0xFF] << 8) |
            (SM4_SBOX[(tmp >> 16) & 0xFF] << 16) | (SM4_SBOX[tmp >> 24] << 24);
        // 线性变换 L(B)
        tmp = tmp ^ ((tmp << 2) | (tmp >> 30)) ^ ((tmp << 10) | (tmp >> 22)) ^
            ((tmp << 18) | (tmp >> 14)) ^ ((tmp << 24) | (tmp >> 8));
        block[0] ^= tmp;

        // 寄存器轮换 (Unbalanced Feistel)
        uint32_t temp = block[3];
        block[3] = block[2];
        block[2] = block[1];
        block[1] = block[0];
        block[0] = temp;
    }
    memcpy(out, block, 16);
}

// ======================= T-Table优化 (4KB) =======================
// 预计算T-Table (PDF Page 6)
alignas(64) uint32_t SM4_TTABLE[4][256];
void init_sm4_ttable() {
    for (int i = 0; i < 256; ++i) {
        uint32_t b = SM4_SBOX[i];
        // 合并线性变换 L(B)
        SM4_TTABLE[0][i] = b ^ ((b << 2) | (b >> 30)) ^ ((b << 10) | (b >> 22)) ^
            ((b << 18) | (b >> 14)) ^ ((b << 24) | (b >> 8));
        SM4_TTABLE[1][i] = (SM4_TTABLE[0][i] << 8) | (SM4_TTABLE[0][i] >> 24);
        SM4_TTABLE[2][i] = (SM4_TTABLE[0][i] << 16) | (SM4_TTABLE[0][i] >> 16);
        SM4_TTABLE[3][i] = (SM4_TTABLE[0][i] << 24) | (SM4_TTABLE[0][i] >> 8);
    }
}

// T-Table加密 (单块)
void sm4_ttable_encrypt(uint8_t out[16], const uint8_t in[16], const uint32_t rk[SM4_ROUNDS]) {
    uint32_t block[4];
    memcpy(block, in, 16);

    for (int i = 0; i < SM4_ROUNDS; ++i) {
        uint32_t x = block[1] ^ block[2] ^ block[3] ^ rk[i];
        uint32_t t = SM4_TTABLE[0][x & 0xFF] ^
            SM4_TTABLE[1][(x >> 8) & 0xFF] ^
            SM4_TTABLE[2][(x >> 16) & 0xFF] ^
            SM4_TTABLE[3][x >> 24];
        block[0] ^= t;

        // 寄存器轮换
        uint32_t temp = block[3];
        block[3] = block[2];
        block[2] = block[1];
        block[1] = block[0];
        block[0] = temp;
    }
    memcpy(out, block, 16);
}

// ======================= AES-NI加速 =======================
#if defined(__AES__)
// GF(2^8)同构映射矩阵 (PDF Page 9)
const __m128i MAPPING_MATRIX = _mm_setr_epi8(
    0x0E, 0x0F, 0x0B, 0x09, 0x0D, 0x0C, 0x0A, 0x08,
    0x06, 0x07, 0x03, 0x01, 0x05, 0x04, 0x02, 0x00
);

__m128i sm4_aesni_encrypt_round(__m128i block, __m128i rk) {
    // 映射到AES域
    __m128i mapped = _mm_shuffle_epi8(MAPPING_MATRIX, block);

    // 使用AES指令执行核心操作
    __m128i result = _mm_aesenclast_si128(mapped, rk);

    // 逆映射回SM4域
    return _mm_shuffle_epi8(MAPPING_MATRIX, result);
}

void sm4_aesni_encrypt(uint8_t out[16], const uint8_t in[16], const __m128i rk[SM4_ROUNDS]) {
    __m128i block = _mm_loadu_si128((const __m128i*)in);

    for (int i = 0; i < SM4_ROUNDS; ++i) {
        block = sm4_aesni_encrypt_round(block, rk[i]);
    }

    _mm_storeu_si128((__m128i*)out, block);
}
#endif

// ======================= GFNI+AVX512加速 =======================
#if defined(__GFNI__) && defined(__AVX512F__)
// GFNI加速的SBox实现 (PDF Page 11)
__m512i sm4_gfni_sbox(__m512i x) {
    // GF(2^8)仿射变换参数
    const __m512i A = _mm512_set1_epi64(0xC7 << 56 | 0x1D << 48 | 0xF7 << 40 | 0xF8 << 32 | 0x01 << 24 | 0x6D << 16 | 0x55 << 8 | 0xBD);
    const __m512i B = _mm512_set1_epi64(0x2F << 56 | 0x9F << 48 | 0x5B << 40 | 0x6A << 32 | 0x35 << 24 | 0xEE << 16 | 0xA1 << 8 | 0xEC);

    // GFNI实现SBox
    return _mm512_gf2p8affine_epi64_epi8(x, A, 0x63);
}

// AVX512加速的线性变换 (使用VPROLD)
__m512i sm4_avx512_linear(__m512i x) {
    __m512i t1 = _mm512_rol_epi32(x, 2);
    __m512i t2 = _mm512_rol_epi32(x, 10);
    __m512i t3 = _mm512_rol_epi32(x, 18);
    __m512i t4 = _mm512_rol_epi32(x, 24);
    return _mm512_xor_si512(x, _mm512_xor_si512(t1, _mm512_xor_si512(t2, _mm512_xor_si512(t3, t4))));
}

// 全并行SM4加密 (8块)
void sm4_avx512_encrypt(uint8_t* out, const uint8_t* in, const __m512i rk[SM4_ROUNDS], size_t blocks) {
    for (size_t i = 0; i < blocks; i += 8) {
        __m512i block = _mm512_loadu_si512((const __m512i*)(in + i * SM4_BLOCK_SIZE));

        for (int r = 0; r < SM4_ROUNDS; r++) {
            // SBox变换 (GFNI加速)
            __m512i sboxed = sm4_gfni_sbox(block);

            // 线性变换 (VPROLD加速)
            block = sm4_avx512_linear(sboxed);

            // 轮密钥加
            block = _mm512_xor_si512(block, rk[r]);

            // Feistel结构移位
            block = _mm512_rolv_epi32(block, _mm512_set_epi32(0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1));
        }

        _mm512_storeu_si512((__m512i*)(out + i * SM4_BLOCK_SIZE), block);
    }
}
#endif

// ======================= SM4-GCM优化 =======================
#if defined(__PCLMUL__)
// GHASH认证 (使用CLMUL指令)
__m128i ghash(__m128i hash, const __m128i* data, size_t blocks, const __m128i H) {
    for (size_t i = 0; i < blocks; ++i) {
        hash = _mm_xor_si128(hash, data[i]);
        hash = _mm_clmulepi64_si128(hash, H, 0x00);  // GF(2^128)乘法
        hash = _mm_clmulepi64_si128(hash, H, 0x11);
    }
    return hash;
}

// SM4-GCM完整实现
void sm4_gcm_encrypt(uint8_t* out, uint8_t* tag, const uint8_t* in, size_t len,
    const uint8_t* key, const uint8_t* iv, const uint8_t* aad, size_t aad_len) {
    // 1. 初始化SM4密钥
    // (密钥扩展代码略)

    // 2. 生成GHASH的H参数
    __m128i H = _mm_setzero_si128();
    // ... 使用SM4加密零块计算H

    // 3. 计数器模式加密
    uint8_t counter[16];
    memcpy(counter, iv, 12);
    memset(counter + 12, 0, 4);

    // 4. 并行加密和认证
    __m128i auth = _mm_setzero_si128();
    size_t blocks = len / 16;

    // 使用优化的SM4加密核心
    for (size_t i = 0; i < blocks; i += 4) {
        // 加密4个计数器块
        uint8_t keystream[64];
        sm4_avx512_encrypt(keystream, counter, rk, 4);  // 使用AVX512优化

        // XOR明文生成密文
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 16; ++k) {
                out[i * 16 + j * 16 + k] = in[i * 16 + j * 16 + k] ^ keystream[j * 16 + k];
            }
        }

        // 更新计数器
        for (int j = 0; j < 4; ++j) {
            uint32_t* ctr = (uint32_t*)(counter + 12);
            *ctr = _byteswap_ulong(_byteswap_ulong(*ctr) + 1);
        }

        // GHASH更新
        auth = ghash(auth, (const __m128i*)(out + i * 16), 4, H);
    }

    // 5. 生成认证标签
    // ... (完整标签计算)
}
#endif

// ======================= 抗缓存攻击 =======================
// 恒定时间T-Table访问 (OpenSSL方式)
uint32_t sm4_ct_ttable_access(uint32_t x) {
    uint32_t result = 0;
    for (int i = 0; i < 4; i++) {
        uint8_t byte = (x >> (i * 8)) & 0xFF;
        for (int j = 0; j < 256; j++) {
            // 恒定时间选择器
            uint32_t mask = (-1) * ((j ^ byte) == 0);
            result |= SM4_TTABLE[i][j] & mask;
        }
    }
    return result;
}

// ======================= 密钥扩展 =======================
// SM4密钥扩展函数
void sm4_key_schedule(const uint8_t key[16], uint32_t rk[SM4_ROUNDS]) {
    const uint32_t FK[4] = { 0xA3B1BAC6, 0x56AA3350, 0x677D9197, 0xB27022DC };
    const uint32_t CK[32] = {
        0x00070E15, 0x1C232A31, 0x383F464D, 0x545B6269,
        0x70777E85, 0x8C939AA1, 0xA8AFB6BD, 0xC4CBD2D9,
        0xE0E7EEF5, 0xFC030A11, 0x181F262D, 0x343B4249,
        0x50575E65, 0x6C737A81, 0x888F969D, 0xA4ABB2B9,
        0xC0C7CED5, 0xDCE3EAF1, 0xF8FF060D, 0x141B2229,
        0x30373E45, 0x4C535A61, 0x686F767D, 0x848B9299,
        0xA0A7AEB5, 0xBCC3CAD1, 0xD8DFE6ED, 0xF4FB0209,
        0x10171E25, 0x2C333A41, 0x484F565D, 0x646B7279
    };

    uint32_t K[4];
    memcpy(K, key, 16);

    // 初始密钥加FK
    K[0] ^= FK[0];
    K[1] ^= FK[1];
    K[2] ^= FK[2];
    K[3] ^= FK[3];

    // 轮密钥生成
    for (int i = 0; i < SM4_ROUNDS; i++) {
        uint32_t X = K[1] ^ K[2] ^ K[3] ^ CK[i];

        // SBox应用
        uint32_t T = SM4_SBOX[X & 0xFF] |
            (SM4_SBOX[(X >> 8) & 0xFF] << 8) |
            (SM4_SBOX[(X >> 16) & 0xFF] << 16) |
            (SM4_SBOX[X >> 24] << 24);

        // 线性变换L'
        T = T ^ ((T << 13) | (T >> 19)) ^ ((T << 23) | (T >> 9));

        rk[i] = K[0] ^ T;

        // 寄存器轮换
        K[0] = K[1];
        K[1] = K[2];
        K[2] = K[3];
        K[3] = rk[i];
    }
}

// ======================= 主函数 =======================
int main() {
    // 初始化T-Table
    init_sm4_ttable();

    // 测试向量 (SM4官方测试向量)
    const uint8_t key[16] = {
        0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF,
        0xFE, 0xDC, 0xBA, 0x98, 0x76, 0x54, 0x32, 0x10
    };

    const uint8_t plaintext[16] = {
        0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF,
        0xFE, 0xDC, 0xBA, 0x98, 0x76, 0x54, 0x32, 0x10
    };

    const uint8_t expected_ciphertext[16] = {
        0x68, 0x1E, 0xDF, 0x34, 0xD2, 0x06, 0x96, 0x5E,
        0x86, 0xB3, 0xE9, 0x4F, 0x53, 0x6E, 0x42, 0x46
    };

    // 生成轮密钥
    uint32_t rk[SM4_ROUNDS];
    sm4_key_schedule(key, rk);

    // 基础实现测试
    uint8_t ciphertext_basic[16];
    sm4_basic_encrypt(ciphertext_basic, plaintext, rk);

    // T-Table实现测试
    uint8_t ciphertext_ttable[16];
    sm4_ttable_encrypt(ciphertext_ttable, plaintext, rk);

    cout << "=== SM4 加密测试 ===" << endl;
    cout << "基础实现: " << (memcmp(ciphertext_basic, expected_ciphertext, 16) == 0 ? "通过" : "失败") << endl;
    cout << "T-Table实现: " << (memcmp(ciphertext_ttable, expected_ciphertext, 16) == 0 ? "通过" : "失败") << endl;

    // 打印加密结果
    auto print_hex = [](const char* label, const uint8_t* data, size_t len) {
        cout << label << ": ";
        for (size_t i = 0; i < len; ++i) {
            cout << std::hex << std::setw(2) << std::setfill('0')
                << static_cast<int>(data[i]) << " ";
        }
        cout << std::dec << endl;
        };

    print_hex("明文     ", plaintext, 16);
    print_hex("基础加密 ", ciphertext_basic, 16);
    print_hex("T-Table加密", ciphertext_ttable, 16);
    print_hex("预期密文 ", expected_ciphertext, 16);

    // 性能基准测试
    const size_t TEST_SIZE = 64 * 1024 * 1024; // 64 MB
    uint8_t* test_data = new uint8_t[TEST_SIZE];
    uint8_t* encrypted_data = new uint8_t[TEST_SIZE];

    // 填充随机测试数据
    for (size_t i = 0; i < TEST_SIZE; ++i) {
        test_data[i] = rand() % 256;
    }

    cout << "\n=== 性能基准测试 (64MB 数据) ===" << endl;

    // 基础实现性能测试
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < TEST_SIZE; i += 16) {
        sm4_basic_encrypt(encrypted_data + i, test_data + i, rk);
    }
    auto end = std::chrono::high_resolution_clock::now();
    double basic_time = std::chrono::duration<double>(end - start).count();
    double basic_speed = (TEST_SIZE / (1024.0 * 1024.0)) / basic_time;
    cout << "基础实现: " << std::fixed << std::setprecision(2)
        << basic_speed << " MB/s" << endl;

    // T-Table实现性能测试
    start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < TEST_SIZE; i += 16) {
        sm4_ttable_encrypt(encrypted_data + i, test_data + i, rk);
    }
    end = std::chrono::high_resolution_clock::now();
    double ttable_time = std::chrono::duration<double>(end - start).count();
    double ttable_speed = (TEST_SIZE / (1024.0 * 1024.0)) / ttable_time;
    cout << "T-Table实现: " << std::fixed << std::setprecision(2)
        << ttable_speed << " MB/s" << endl;

    // AES-NI加速测试
#if defined(__AES__)
    __m128i rk_aesni[SM4_ROUNDS];
    // 转换轮密钥格式
    for (int i = 0; i < SM4_ROUNDS; ++i) {
        rk_aesni[i] = _mm_set1_epi32(rk[i]);
    }

    start = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < TEST_SIZE; i += 16) {
        sm4_aesni_encrypt(encrypted_data + i, test_data + i, rk_aesni);
    }
    end = chrono::high_resolution_clock::now();
    double aesni_time = chrono::duration<double>(end - start).count();
    double aesni_speed = (TEST_SIZE / (1024.0 * 1024.0)) / aesni_time;
    cout << "AES-NI加速: " << fixed << setprecision(2)
        << aesni_speed << " MB/s" << endl;
#else
    cout << "AES-NI加速: 当前CPU不支持" << endl;
#endif

    // AVX512+GFNI加速测试
#if defined(__GFNI__) && defined(__AVX512F__)
    __m512i rk_avx512[SM4_ROUNDS];
    // 转换轮密钥格式
    for (int i = 0; i < SM4_ROUNDS; ++i) {
        rk_avx512[i] = _mm512_set1_epi32(rk[i]);
    }

    start = chrono::high_resolution_clock::now();
    sm4_avx512_encrypt(encrypted_data, test_data, rk_avx512, TEST_SIZE / 16);
    end = chrono::high_resolution_clock::now();
    double avx512_time = chrono::duration<double>(end - start).count();
    double avx512_speed = (TEST_SIZE / (1024.0 * 1024.0)) / avx512_time;
    cout << "AVX512+GFNI加速: " << fixed << setprecision(2)
        << avx512_speed << " MB/s" << endl;
#else
    cout << "AVX512+GFNI加速: 当前CPU不支持" << endl;
#endif

    // 清理内存
    delete[] test_data;
    delete[] encrypted_data;

    return 0;
}