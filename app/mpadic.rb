module X
  class Mpadic

    ############################## CLASS MEMBERS ##############################
    attr_accessor  :primes, :primes_exp, :r, :padics

    ############################## CONSTRUCTOR ##############################
    # c = numerator
    # d = denominator
    def initialize(primes,r,c=0,d=1)
      @primes = primes
      @r = r

      @padics = []
      @primes_exp = []

      primes.each do |p|
        @padics << Padic.new(p,r,c,d)
        @primes_exp << p**r
      end
    end

    ############################## OBJECT METHODS ##############################
    def decode
      # Compute p
      p = 1
      (0..primes.size-1).each do |i|
        p = p * primes_exp[i]
      end

      # Compute residue
      i = 0
      residue = 0
      padics.each do |padic|
        p_factor = p / (primes_exp[i])
        p_prime = Tools.mod_inverse(p_factor,primes_exp[i])
        residue = residue + ((p_factor * p_prime * padic.hensel_code))
        i+=1
      end

      residue = residue % p

      Tools.ctr(p, residue , Math.sqrt(p))
    end

    def add(b)
      result = []
      (0..self.padics.size-1).each do |i|
        result << self.padics[i].add(b.padics[i])
      end

      c = Mpadic.new(self.primes, self.r)
      c.padics = result
      c
    end

    # def add_parallel(b)
    #   result = []
    #   # Parallel.each(0..self.padics.size-1, :in_threads => 2) { |i| result[i] = self.padics[i].add(b.padics[i]) }
    #   c = Mpadic.new(self.primes, self.r)
    #   c.padics = result
    #   c
    # end

    def sub(b)
      result = []
      (0..self.padics.size-1).each do |i|
        result << self.padics[i].sub(b.padics[i])
      end

      c = Mpadic.new(self.primes, self.r)
      c.padics = result
      c
    end

    def mul(b)
      result = []
      (0..self.padics.size-1).each do |i|
        result << self.padics[i].mul(b.padics[i])
      end

      c = Mpadic.new(self.primes, self.r)
      c.padics = result
      c
    end

    def mul_parallel(b)
      result = []

      r = Parallel.each(0..self.padics.size-1, :in_threads => 4) { |i| result[i] = self.padics[i].mul(b.padics[i]) }

      c = Mpadic.new(self.primes, self.r)
      c.padics = result
      c
    end

    def div(b)
      result = []
      (0..self.padics.size-1).each do |i|
        result << self.padics[i].div(b.padics[i])
      end

      c = Mpadic.new(self.primes, self.r)
      c.padics = result
      c
    end

    def negate
      self.mul(MPadic.new(primes,r,-1,1))
    end


    # Prints value
    def to_s
      padics
    end

    def inspect
      to_s
    end


  end

end
