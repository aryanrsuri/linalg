pub const std = @import("std");
pub const rgen = std.rand.DefaultPrng;
pub var rand = rgen.init(0);
pub fn Matrix(comptime T: type) type {
    return struct {
        const Row = struct { cols: usize, elements: *T };
        const Self = @This();
        rows: usize,
        cols: usize,
        elements: []T,
        region: std.mem.Allocator,

        pub fn alloc(allocator: std.mem.Allocator, comptime rows: usize, comptime cols: usize) Self {
            const elements = allocator.alloc(T, rows * cols) catch @panic("Allocation failed");
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

        pub fn rmask(self: *Self, radix: T) void {
            var i: usize = 0;
            while (i < self.rows) : (i += 1) {
                var j: usize = 0;
                while (j < self.cols) : (j += 1) {
                    if (@typeInfo(T) == .Float) {
                        self.set(i, j, rand.random().float(T) * radix);
                    } else {
                        self.set(i, j, rand.random().int(T) * radix);
                    }
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

        /// IN PROGRESS DO NOT USE
        pub fn submatrix(self: *Self, i: usize, j: usize, M: *Matrix(T)) !void {
            var iter: usize = i;
            var jter: usize = j;
            while (iter < j) : (iter += 1) {
                std.debug.print("{} {}\n", .{ iter, jter });
                const index: usize = @as(usize, self.elements[iter]);
                M.elements[index];
                jter += 1;
            }
        }

        pub fn dot(self: *Self, A: *Matrix(T), M: *Matrix(T)) !void {
            if (self.rows != A.rows) return error.MatrixSpaceError;
            if (A.rows != M.rows) return error.MatrixSpaceError;
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

        pub fn trace(self: *Self) !T {
            if (self.rows != self.cols) return error.NotSquare;
            var i: usize = 0;
            var tr: T = undefined;
            while (i < self.rows) : (i += 1) {
                var j: usize = 0;
                while (j < self.cols) : (j += 1) {
                    if (i == j) tr += self.get(i, j);
                }
            }

            return tr;
        }
        pub fn determinant() T {}
        pub fn magnitude() T {}
        pub fn inverse() T {}

        pub fn transpose(self: *Self, M: *Matrix(T)) !void {
            if (self.rows != M.cols) return error.MatrixSpaceError;
            var i: usize = 0;
            while (i < self.rows) : (i += 1) {
                var j: usize = 0;
                while (j < self.cols) : (j += 1) {
                    var v: T = self.get(i, j);
                    M.set(j, i, v);
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

        pub fn copy(self: *Self, M: *Matrix(T)) void {
            if (self.rows != A.rows) return error.MatrixSpaceError;
            var i: usize = 0;
            while (i < self.rows) : (i += 1) {
                var j: usize = 0;
                while (j < self.cols) : (j += 1) {
                    const v = self.get(i, j);
                    M.set(i, j, v);
                }
            }
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
            var i: usize = 0;
            while (i < self.rows) : (i += 1) {
                var j: usize = 0;
                while (j < self.cols) : (j += 1) {
                    if (i == j) {
                        self.set(i, j, 1);
                    } else self.set(i, j, 0);
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
    var I = Matrix(i64).alloc(T, 3, 2);
    var IrE = Matrix(i64).alloc(T, 2, 3);
    var A = Matrix(f64).alloc(T, 2, 2);
    var M = Matrix(f64).alloc(T, 2, 2);
    var E = Matrix(f64).alloc(T, 2, 3);
    var TrE = Matrix(f64).alloc(T, 3, 2);

    E.rmask(1);
    I.rmask(1);
    E.print("E");
    I.print("I");
    try E.transpose(&TrE);
    try I.transpose(&IrE);

    TrE.print("E_tr");
    IrE.print("I_tr");
    var two = Mat2(f64, T);
    var TH = Mat3(f64, T);
    TH.ones();
    TH.print("m3");
    TH.set(1, 1, 100.0);
    TH.print("m3");
    TH.rmask(10);
    TH.print("m3");
    defer {
        R.dealloc();
        A.dealloc();
        I.dealloc();
        M.dealloc();
        E.dealloc();
        TrE.dealloc();
        IrE.dealloc();
        TH.dealloc();
        two.dealloc();
    }

    M.mask(3.0);
    std.debug.print("\n\t Trace of Mat3rand {any}\n", .{TH.trace()});
    R.mask(5.0);
    try M.dot(&R, &A);
    R.scale(2.0);

    R.print("R");
    // try R.submatrix(0, 3, &two);
    // two.print("2x2 subof R");
    A.print("A");
    M.print("M");
    M.identity();
    M.print("I");
}
