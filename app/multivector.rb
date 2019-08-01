module X
  class Multivector

    attr_accessor :obj_type

    def plaintext?
      obj_type == "plaintext"
    end

    def ciphertext?
      obj_type == "ciphertext"
    end

    def dump
      Marshal.dump(self)
    end

    def save_to_disk
      filename = ("%10.9f" % Time.now).to_s
      file = open(Dir.pwd + "/saved/" + filename,"wb") do |f|
        f.write dump
      end
      filename
    end

  end
end
