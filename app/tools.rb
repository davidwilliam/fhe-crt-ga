module X
  class Tools

    def self.experimental_hensel_packing(m, p, q, rs_data, delta_data)
      p_bits = p.bit_length

      rs = rs_data

      d1 = Rational(rs[0],rs[1])
      d2 = Rational(rs[2],rs[3])

      data = delta_data

      data << (m - data.inject(:+))

      mm = Multivector2D.new data

      data_ = [
        data[0] + d1,
        data[1] - d1,
        data[2] + d2,
        data[3] - d2
      ]

      mm_ = Multivector2D.new data_

      data_p = data_.map{|d|
        FHE.hensel_encoding(d, p)
      }

      mm_p = Multivector2Dff.new data_p, q

      mm_p
    end

    def self.scan_hensel_packing(m, p, q, rs_data, delta_data, sample)
      total = 0

      rs_data.each do |rs|
        delta_data.each do |delta|
          mm = experimental_hensel_packing(m, p, q, rs, delta)
          if mm.data == sample
            total += 1
          end
        end
      end

      puts "total = #{total}"

      nil
    end

    def self.mod_inverse(num, mod)
      g, a, b = extended_gcd(num, mod)
      unless g == 1
        raise ZeroDivisionError.new("#{num} has no inverse modulo #{mod}")
      end
      a % mod
    end

    def self.extended_gcd(x, y)
      if x < 0
        g, a, b = extended_gcd(-x, y)
        return [g, -a, b]
      end
      if y < 0
        g, a, b = extended_gcd(x, -y)
        return [g, a, -b]
      end
      r0, r1 = x, y
      a0 = b1 = 1
      a1 = b0 = 0
      until r1.zero?
        q = r0 / r1
        r0, r1 = r1, r0 - q*r1
        a0, a1 = a1, a0 - q*a1
        b0, b1 = b1, b0 - q*b1
      end
      [r0, a0, b0]
    end

    def self.ctr (a0, a1, n)
      a = []
      a << a0
      a << a1

      y = []
      y << 0
      y << 1

      i = 1
      while(a[i] > n)
        q = (a[i-1]/a[i]).floor
        a << a[i-1] - q*a[i]
        y << y[i-1] + q*y[i]
        i = i + 1
      end

      Rational((-1)**(i+1) * a[i],y[i])
    end

  end

end
