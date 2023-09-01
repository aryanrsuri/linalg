const std = @import("std");
const rgen = std.rand.DefaultPrng;

pub fn Matrix(comptime T: type) type {
    return struct {
        const Row = struct { cols: usize, elements: *T };
        const Self = @This();
        rows: usize,
        cols: usize,
        elements: []T,
        region: std.mem.Allocator,

        pub fn alloc(allocator: std.mem.Allocator, comptime rows: usize, comptime cols: usize) Self {
            const elements = allocator.alloc(T, rows * cols) catch @panic("Allocaiton failed");
            return .{
                .rows = rows,
                .cols = cols,
                .region = allocator,
                .elements = elements,
            };
        }

        pub fn dealloc(self: *Self) void {
            self.region.free(self.elements);
            self.* = undefined;
        }

        pub fn set(self: *Self, i: usize, j: usize, scalar: T) void {
            self.elements[(i) * (self.cols) + (j)] = scalar;
        }

        pub fn get(self: *Self, i: usize, j: usize) T {
            return self.elements[(i) * (self.cols) + (j)];
        }

        pub fn mask(self: *Self, scalar: T) void {
            var i: usize = 0;
            while (i < self.rows) : (i += 1) {
                var j: usize = 0;
                while (j < self.cols) : (j += 1) {
                    self.set(i, j, scalar);
                }
            }
        }

        pub fn rmask(self: *Self) void {
            var rand = rgen.init(0);
            var i: usize = 0;
            while (i < self.rows) : (i += 1) {
                var j: usize = 0;
                while (j < self.cols) : (j += 1) {
                    // if (tinfo == std.builtin.Type.Int) {
                    // self.set(i, j, rand.random().int(T));
                    // } else {
                    self.set(i, j, rand.random().float(T));
                    // }
                }
            }
        }
        pub fn scale(self: *Self, scalar: T) void {
            var i: usize = 0;
            while (i < self.rows) : (i += 1) {
                var j: usize = 0;
                while (j < self.cols) : (j += 1) {
                    const v = self.get(i, j);
                    self.set(i, j, scalar * v);
                }
            }
        }

        pub fn dot(self: *Self, A: *Matrix(T), M: *Matrix(T)) void {
            var i: usize = 0;
            while (i < self.rows) : (i += 1) {
                var j: usize = 0;
                while (j < self.cols) : (j += 1) {
                    var k: usize = 0;
                    var v: T = 0;
                    M.set(i, j, 0);
                    while (k < self.cols) : (k += 1) {
                        v += self.get(i, k) * A.get(i, k);
                        M.set(i, j, v);
                    }
                }
            }
        }

        pub fn trace() T {}
        pub fn determinant() T {}
        pub fn inverse() T {}

        pub fn transpose(self: *Self, M: *Matrix(T)) void {
            var i: usize = 0;
            while (i < self.rows) : (i += 1) {
                var j: usize = 0;
                while (j < self.cols) : (j += 1) {
                    var v: T = self.get(i, j);
                    M.set(i, j, v);
                }
            }
        }

        pub fn sum(self: *Self, A: *Matrix(T), M: *Matrix(T)) void {
            var i: usize = 0;
            while (i < self.rows) : (i += 1) {
                var j: usize = 0;
                while (j < self.cols) : (j += 1) {
                    var k: usize = 0;
                    var v: T = 0;
                    M.set(i, j, 0);
                    while (k < self.cols) : (k += 1) {
                        v += self.get(i, k) + A.get(i, k);
                        M.set(i, j, v);
                    }
                }
            }
        }

        pub fn copy(self: *Self, A: *Matrix(T)) void {
            var i: usize = 0;
            while (i < self.rows) : (i += 1) {
                var j: usize = 0;
                while (j < self.cols) : (j += 1) {
                    const v = self.get(i, j);
                    A.set(i, j, v);
                }
            }
        }

        pub fn print(self: *Self, title: []const u8) void {
            var i: usize = 0;
            std.debug.print("\n{s} = [\n", .{title});
            while (i < self.rows) : (i += 1) {
                var j: usize = 0;
                std.debug.print("\t", .{});
                while (j < self.cols) : (j += 1) {
                    std.debug.print("{any} ", .{self.get(i, j)});
                }
                std.debug.print("\n", .{});
            }
            std.debug.print(" ] \n ", .{});
        }
        // Alias functions //

        pub fn zeros(self: *Self) void {
            self.mask(0);
        }

        pub fn ones(self: *Self) void {
            self.mask(1);
        }

        pub fn identity(self: *Self) void {
            if (self.rows != self.cols) @panic("Matrix not square");
        }
    };
}

pub fn Mat2(comptime T: type, alloc: std.mem.Allocator) Matrix(T) {
    return Matrix(T).alloc(alloc, 2, 2);
}
pub fn Mat3(comptime T: type, alloc: std.mem.Allocator) Matrix(T) {
    return Matrix(T).alloc(alloc, 3, 3);
}
pub fn Mat4(comptime T: type, alloc: std.mem.Allocator) Matrix(T) {
    return Matrix(T).alloc(alloc, 3, 3);
}

test "Matrix" {
    const T = std.testing.allocator;
    var R = Matrix(f64).alloc(T, 2, 2);
    var A = Matrix(f64).alloc(T, 2, 2);
    var M = Matrix(f64).alloc(T, 2, 2);
    var E = Matrix(f64).alloc(T, 2, 3);
    var TrE = Matrix(f64).alloc(T, 3, 2);
    E.rmask();
    E.print("E");
    E.transpose(&TrE);

    TrE.print("E");
    var TH = Mat3(f64, T);
    TH.ones();
    std.debug.print("R\n {any}\n", .{TH.elements});
    TH.set(1, 1, 100.0);
    std.debug.print("R\n {any}\n", .{TH.elements});
    defer {
        R.dealloc();
        A.dealloc();
        M.dealloc();
        E.dealloc();
        TrE.dealloc();
        TH.dealloc();
    }

    M.mask(3.0);
    R.mask(5.0);
    M.dot(&R, &A);
    R.scale(2.0);

    R.print("R");
    A.print("A");
    M.print("M");
    M.print("M_t");
}
