module X
  class Padic

    ############################## CLASS MEMBERS ##############################
    attr_accessor  :p, :r, :digits


    ############################## CONSTRUCTOR ##############################
    # c = numberator
    # d = denominator
    def initialize(p,r,c=0,d=1)
      @p = p
      @r = r
      if(c!=0)
        @digits = Padic.digits(p,r,Padic.hensel_code(p,r,c,d))
      else
        @digits = Array.new(r,0)
      end
    end

    ############################## STATIC METHODS ##############################
    # PARAMETERS: prime p, digit precision r, rational a=c/d
    # RETURNS:    hensel code, which enables the generation of finite p-adic digits
    def self.hensel_code(p,r,c,d)
      pr = p**r
      gcd,x,y = Tools.extended_gcd(pr,d)
      return m = (c * y) % pr
    end

    # PARAMETERS: prime p, digit precision r, hensel code
    # RETURNS:    finite p-adic digits
    def self.digits(p,r,hensel_code)
      digits = []
      (r-1).times do
        div,mod = hensel_code.divmod(p)
        hensel_code = div
        digits << mod
      end
      digits << hensel_code
      digits
    end

    # PARAMETERS: prime p, digit precision r, p-adic digits
    # RETURNS:    hensel code
    def self.digits_to_hensel_code(p,r,digits)
      i = 0
      hensel_code = 0
      r.times do
        hensel_code += (digits[i] * p**i)
        i += 1
      end
      hensel_code
    end

    ############################## OBJECT METHODS ##############################
    # RETURNS:    hensel code of the p-adic numbers
    def hensel_code
      Padic.digits_to_hensel_code(p,r,digits)
    end

    # RETURNS:    the rational representation of the p-adic number
    # NOTE:       if rational a=c/d, 2cd < p^r
    def decode
      Tools.ctr(p ** r, Padic.digits_to_hensel_code(p, r, digits), Math.sqrt(p**r/2))
    end

    # RETURNS:    self + b
    def add b
      result = []

      carry = 0
      (0..r-1).each do |i|
        carry, mod = (self.digits[i] + b.digits[i] + carry).divmod(self.p)
        result << mod
      end

      c = Padic.new(self.p, self.r)
      c.digits = result
      c

      # result = (self.hensel_code + b.hensel_code) % (self.p**self.r)
      # c = Padic.new(self.p, self.r)
      # c.digits = Padic.digits(self.p, self.r, result)
      # c

    end

    # RETURNS:    self - b
    def sub b
      self.add(b.negate)
    end

    # RETURNS:    self * b
    def mul b
      x = self
      y = b
      result = [0]
      j = 0
      y.digits.each do |m|
        c = 0
        i = j
        x.digits.each do |d|
          v = result[i]
          result << 0 if v.zero?
          c, v = (v + c + d*m).divmod(x.p)
          result[i] = v
          i += 1
        end
        if result[i] != nil
          result[i] += c
        else
          result[i] = 0
        end

        j += 1
      end

      c = Padic.new(self.p, self.r)
      c.digits = result.take(self.r)
      c
    end

    # RETURNS:    self / b
    def div b
      self.mul(b.inverse)
    end

    # RETURNS:    self * -1
    def negate
      self.mul(Padic.new(p,r,-1,1))
    end

    # RETURNS:    self ^ -1
    def inverse
      inv_hensel_code = Tools.mod_inverse(self.hensel_code,p**r)
      c = Padic.new(self.p, self.r)
      c.digits = Padic.digits(self.p, self.r,inv_hensel_code)
      c
    end


    # Prints value
    def to_s
      digits.to_s
    end

    def inspect
      to_s
    end

  end

end
