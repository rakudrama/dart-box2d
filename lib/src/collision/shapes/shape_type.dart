// Copyright 2012 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

/**
 * An enum class with the different kinds of shapes. The two kinds of valid
 * shapes are circles and polygons.
 */

part of box2d;

class ShapeType {
  // Circle and Polygon.
  static const int TYPE_COUNT = 2;

  static const int UNKNOWN = -1;
  static const int CIRCLE = 0;
  static const int POLYGON = 1;
}
