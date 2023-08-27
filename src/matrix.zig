const std = @import("std");

pub fn Matrix(comptime T: type, comptime m: usize, comptime n: usize) T {
    return struct { m, n };
}
