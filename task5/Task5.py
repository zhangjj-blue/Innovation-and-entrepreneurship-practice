import random
from typing import Tuple
from hashlib import sha256

# ===================== SM2 椭圆曲线参数 =====================
# SM2参数 (256-bit prime field)
p_sm2 = 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF
a_sm2 = 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC
b_sm2 = 0x28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93
n_sm2 = 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123
Gx_sm2 = 0x32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7
Gy_sm2 = 0xBC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0

# ===================== 比特币 secp256k1 曲线参数 =====================
p_btc = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
a_btc = 0
b_btc = 7
n_btc = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
Gx_btc = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
Gy_btc = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8


class SM2:
    def __init__(self):
        self.p = p_sm2
        self.a = a_sm2
        self.b = b_sm2
        self.n = n_sm2
        self.G = (Gx_sm2, Gy_sm2)

    # 蒙哥马利模约 (SM2优化版)
    def mont_reduce(self, c: int) -> int:
        # 获取低256位和高256位
        c_low = c & ((1 << 256) - 1)
        c_high = c >> 256

        # 利用SM2素数特性快速约简 (p = 2²⁵⁶ - 2²²⁴ - 2⁹⁶ + 2⁶⁴ - 1)
        # 步骤1: 提取特定比特段
        s1 = (c_low >> 224) << 96  # 取高32位，左移96位
        s2 = (c_low >> 96) & ((1 << 128) - 1)  # 取中间128位
        s3 = (c_low & ((1 << 96) - 1)) << 160  # 取低96位，左移160位

        # 步骤2: 组合计算结果
        t = c_high + s1 - s2 + s3

        # 步骤3: 模约简处理
        if t >= self.p:
            t -= self.p
        if t < 0:
            t += self.p

        return t

    # 点加倍 (雅可比坐标)
    def point_double(self, P: Tuple[int, int, int]) -> Tuple[int, int, int]:
        X1, Y1, Z1 = P
        if Z1 == 0:
            return (0, 0, 0)

        # 优化公式: 2M + 4S + 7A
        T1 = 3 * X1 ** 2 + self.a * Z1 ** 4
        Z3 = 2 * Y1 * Z1
        T2 = 4 * X1 * Y1 ** 2
        X3 = T1 ** 2 - 2 * T2
        Y3 = T1 * (T2 - X3) - 8 * Y1 ** 4
        return (X3 % self.p, Y3 % self.p, Z3 % self.p)

    # 点加 (雅可比坐标)
    def point_add(self, P: Tuple[int, int, int], Q: Tuple[int, int, int]) -> Tuple[int, int, int]:
        X1, Y1, Z1 = P
        X2, Y2, Z2 = Q
        if Z1 == 0:
            return Q
        if Z2 == 0:
            return P

        # 优化公式: 8M + 3S + 7A
        Z1_sq = Z1 ** 2
        Z2_sq = Z2 ** 2
        U1 = X1 * Z2_sq
        U2 = X2 * Z1_sq
        S1 = Y1 * Z2_sq * Z2
        S2 = Y2 * Z1_sq * Z1

        if U1 == U2:
            if S1 != S2:
                return (0, 0, 1)  # 无穷远点
            return self.point_double(P)

        H = U2 - U1
        R = S2 - S1
        H_sq = H ** 2
        H_cu = H * H_sq

        X3 = R ** 2 - H_cu - 2 * U1 * H_sq
        Y3 = R * (U1 * H_sq - X3) - S1 * H_cu
        Z3 = H * Z1 * Z2
        return (X3 % self.p, Y3 % self.p, Z3 % self.p)

    # NAF编码 (非相邻形式)
    def naf_encode(self, k: int, w: int = 5) -> list:
        naf = []
        while k > 0:
            if k & 1:
                ki = 2 - (k % (1 << w))
                k -= ki
            else:
                ki = 0
            naf.append(ki)
            k //= 2
        return naf[::-1]

    # 标量乘法 (NAF优化)
    def scalar_mult(self, k: int, P: Tuple[int, int]) -> Tuple[int, int]:
        # 转换为雅可比坐标
        P_jac = (P[0], P[1], 1)

        # 预计算表 (w=5)
        precomputed = {}
        current = P_jac
        for i in range(-15, 16, 2):
            if i != 0:
                precomputed[i] = current
            current = self.point_add(current, P_jac)

        # NAF编码
        naf = self.naf_encode(k)
        R = (0, 0, 0)  # 无穷远点

        # 滑动窗口乘法
        for bit in naf:
            R = self.point_double(R)
            if bit != 0:
                R = self.point_add(R, precomputed[bit])

        # 转回仿射坐标
        if R[2] == 0:
            return (0, 0)
        z_inv = pow(R[2], self.p - 2, self.p)
        x = R[0] * z_inv ** 2 % self.p
        y = R[1] * z_inv ** 3 % self.p
        return (x, y)

    # SM2签名
    def sign(self, dA: int, msg: bytes) -> tuple[int, int] | None:
        e = int.from_bytes(msg, 'big') % self.n
        while True:
            k = random.randint(1, self.n - 1)
            # 标量乘法返回仿射坐标，我们只需要x坐标
            x1, _ = self.scalar_mult(k, self.G)
            r = (e + x1) % self.n
            if r == 0 or r + k == self.n:
                continue
            s = (pow(1 + dA, self.n - 2, self.n) * (k - r * dA)) % self.n
            if s != 0:
                return (r, s)

    # SM2验签
    def verify(self, PA: Tuple[int, int], msg: bytes, sig: Tuple[int, int]) -> bool:
        r, s = sig
        if not (1 <= r < self.n and 1 <= s < self.n):
            return False

        e = int.from_bytes(msg, 'big') % self.n
        t = (r + s) % self.n
        if t == 0:
            return False

        # 双点乘优化
        sG = self.scalar_mult(s, self.G)  # 返回仿射坐标
        tPA = self.scalar_mult(t, PA)  # 返回仿射坐标

        # 将仿射坐标转换为雅可比坐标
        sG_jac = (sG[0], sG[1], 1)
        tPA_jac = (tPA[0], tPA[1], 1)

        # 点加操作
        result_jac = self.point_add(sG_jac, tPA_jac)

        # 转回仿射坐标
        if result_jac[2] == 0:
            x1 = 0
        else:
            z_inv = pow(result_jac[2], self.p - 2, self.p)
            x1 = result_jac[0] * z_inv ** 2 % self.p

        R = (e + x1) % self.n
        return R == r


class ECDSA:
    """比特币使用的ECDSA实现 (secp256k1曲线)"""

    def __init__(self):
        self.p = p_btc
        self.a = a_btc
        self.b = b_btc
        self.n = n_btc
        self.G = (Gx_btc, Gy_btc)

    # 点加倍 (仿射坐标)
    def point_double(self, P: Tuple[int, int]) -> Tuple[int, int]:
        if P[1] == 0:  # 无穷远点
            return (0, 0)
        x, y = P
        s = (3 * x * x + self.a) * pow(2 * y, self.p - 2, self.p) % self.p
        x3 = (s * s - 2 * x) % self.p
        y3 = (s * (x - x3) - y) % self.p
        return (x3, y3)

    # 点加 (仿射坐标)
    def point_add(self, P: Tuple[int, int], Q: Tuple[int, int]) -> Tuple[int, int]:
        if P[1] == 0:  # P是无穷远点
            return Q
        if Q[1] == 0:  # Q是无穷远点
            return P
        if P[0] == Q[0]:
            if P[1] == Q[1]:
                return self.point_double(P)
            else:  # 点互为逆元
                return (0, 0)
        x1, y1 = P
        x2, y2 = Q
        s = (y2 - y1) * pow(x2 - x1, self.p - 2, self.p) % self.p
        x3 = (s * s - x1 - x2) % self.p
        y3 = (s * (x1 - x3) - y1) % self.p
        return (x3, y3)

    # 标量乘法 (双倍-加法)
    def scalar_mult(self, k: int, P: Tuple[int, int]) -> Tuple[int, int]:
        # 处理无穷远点
        if k % self.n == 0:
            return (0, 0)

        # 将k转为二进制
        bits = bin(k)[2:]
        R = (0, 0)  # 初始化为无穷远点

        for bit in bits:
            R = self.point_double(R)
            if bit == '1':
                R = self.point_add(R, P)
        return R

    # ECDSA签名
    def sign(self, d: int, msg: bytes, k: int = None) -> Tuple[int, int]:
        e = int.from_bytes(sha256(msg).digest(), 'big') % self.n
        k = k or random.randint(1, self.n - 1)
        x, _ = self.scalar_mult(k, self.G)
        r = x % self.n
        s = (pow(k, self.n - 2, self.n) * (e + d * r)) % self.n
        return (r, s)

    # 验证签名
    def verify(self, Q: Tuple[int, int], msg: bytes, sig: Tuple[int, int]) -> bool:
        r, s = sig
        if not (1 <= r < self.n and 1 <= s < self.n):
            return False

        e = int.from_bytes(sha256(msg).digest(), 'big') % self.n
        w = pow(s, self.n - 2, self.n)
        u1 = (e * w) % self.n
        u2 = (r * w) % self.n

        # 计算 u1*G + u2*Q
        P1 = self.scalar_mult(u1, self.G)
        P2 = self.scalar_mult(u2, Q)
        x, _ = self.point_add(P1, P2)

        return r == x % self.n

    # 伪造签名 (使用重复k值)
    def forge_signature(self, sig1: Tuple[int, int], msg1: bytes,
                        sig2: Tuple[int, int], msg2: bytes) -> int:
        r1, s1 = sig1
        r2, s2 = sig2
        if r1 != r2:
            raise ValueError("不同的r值，无法利用漏洞")

        e1 = int.from_bytes(sha256(msg1).digest(), 'big') % self.n
        e2 = int.from_bytes(sha256(msg2).digest(), 'big') % self.n

        # 计算私钥 d = (s2*e1 - s1*e2) / (r*(s1 - s2)) mod n
        numerator = (s2 * e1 - s1 * e2) % self.n
        denominator = (r1 * (s1 - s2)) % self.n
        d = numerator * pow(denominator, self.n - 2, self.n) % self.n
        return d


def sm2_signature_misuse_poc():
    """演示SM2签名误用漏洞"""
    global s_wrong
    print("\n" + "=" * 60)
    print("SM2签名误用漏洞演示")
    print("=" * 60)

    sm2 = SM2()
    dA = 0x1234ABCD  # 测试私钥
    PA = sm2.scalar_mult(dA, sm2.G)
    msg = b"Critical security message"

    # 正确签名
    r_correct, s_correct = sm2.sign(dA, msg)

    # 误用签名 (r计算使用mod p)
    e = int.from_bytes(msg, 'big') % sm2.n
    while True:
        k = random.randint(1, sm2.n - 1)
        x1, _ = sm2.scalar_mult(k, sm2.G)
        r_wrong = (e + x1) % sm2.p  # 误用mod p
        if r_wrong == 0:
            continue
        s_wrong = (pow(1 + dA, sm2.n - 2, sm2.n) * (k - r_wrong * dA)) % sm2.n
        if s_wrong != 0:
            break

    # 验证结果
    valid_correct = sm2.verify(PA, msg, (r_correct, s_correct))
    valid_wrong = sm2.verify(PA, msg, (r_wrong, s_wrong))

    print(f"[正确签名] 验证结果: {valid_correct} (期望True)")
    print(f"[误用签名] 验证结果: {valid_wrong} (期望False)")
    print(f"说明：误用mod p会导致签名验证失败，这是SM2签名算法的正确行为")


def ecdsa_forge_signature_poc():
    """演示ECDSA重复k值攻击"""
    print("\n" + "=" * 60)
    print("ECDSA重复k值攻击演示")
    print("=" * 60)

    ecdsa = ECDSA()

    # 假设中本聪私钥 (未知)
    d_satoshi = 0x1E240  # 测试用

    # 用相同k生成两个签名
    k = 0xABCD1234  # 重复使用的k值
    msg1 = b"Send 10 BTC to Alice"
    msg2 = b"Send 20 BTC to Bob"
    sig1 = ecdsa.sign(d_satoshi, msg1, k)
    sig2 = ecdsa.sign(d_satoshi, msg2, k)

    print(f"消息1: {msg1.decode()}, 签名1: (r={hex(sig1[0])}, s={hex(sig1[1])})")
    print(f"消息2: {msg2.decode()}, 签名2: (r={hex(sig2[0])}, s={hex(sig2[1])})")

    # 伪造私钥
    try:
        d_forged = ecdsa.forge_signature(sig1, msg1, sig2, msg2)
    except ValueError as e:
        print(f"错误: {e}")
        return

    # 验证伪造
    target_msg = b"Transfer all BTC to Attacker"
    forged_sig = ecdsa.sign(d_forged, target_msg)
    valid = ecdsa.verify(ecdsa.scalar_mult(d_forged, ecdsa.G), target_msg, forged_sig)

    print(f"\n真实私钥: {hex(d_satoshi)}")
    print(f"伪造私钥: {hex(d_forged)}")
    print(f"私钥匹配: {d_satoshi == d_forged}")
    print(f"使用伪造私钥签名验证: {'成功' if valid else '失败'}")
    print(f"说明：通过两个使用相同k值的签名，可以计算出私钥")


def main():
    """主函数"""
    print("=" * 60)
    print("椭圆曲线密码学安全演示")
    print("=" * 60)
    print("1. SM2签名误用漏洞演示")
    print("2. ECDSA重复k值攻击演示")
    print("3. 退出")

    while True:
        choice = input("\n请选择演示(1-3): ")

        if choice == '1':
            sm2_signature_misuse_poc()
        elif choice == '2':
            ecdsa_forge_signature_poc()
        elif choice == '3':
            print("程序退出")
            break
        else:
            print("无效选择，请重新输入")


if __name__ == "__main__":
    main()

